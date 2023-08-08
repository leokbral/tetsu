import os
import subprocess
import glob
import time

def running_proteins_structure (_dir, SCREEN, DATABASE, HEADN, ANNOT, num_processors):
    start_alltask_time = time.time()
    def get_num_lines(file_path):
        with open(file_path, "r") as file:
            return sum(1 for _ in file)
    print('*'*64, 'inside running proteins structure')
    os.chdir(_dir+'/content/input')
    if os.path.isfile(f'{_dir}/content/result/result.tab'):
        os.remove(f'{_dir}/content/result/result.tab')
    if os.path.isfile(f'{_dir}/content/result/tmalign_formatted.tab'):
        os.remove(f'{_dir}/content/result/tmalign_formatted.tab')
    if os.path.isfile(f'{_dir}/content/result/fatcat_formatted.tab'):
        os.remove(f'{_dir}/content/result/fatcat_formatted.tab')
    if os.path.isfile(f'{_dir}/content/result/lovoalign_formatted.tab'):
        os.remove(f'{_dir}/content/result/lovoalign_formatted.tab')

   
    if SCREEN == "foldseek":
        for f in os.listdir():
                os.system(f"{_dir}/content/programs/foldseek/bin/foldseek easy-search {f} {_dir}/content/foldseek_data/fs_{DATABASE} {_dir}/content/result/screening/tmp.tab.fmt {_dir}/content/tmpFolder --max-seqs {HEADN} -e inf")
                os.system((f"cut -f 1,2 {_dir}/content/result/screening/tmp.tab.fmt | sort | uniq | perl -ne '@a = split(/\\./, $_); print join(\".\", @a[0 .. $#a-1]).\"\\n\";' > {_dir}/content/result/screening/{f}.tab.fmt"))
                os.remove(f"{_dir}/content/result/screening/tmp.tab.fmt")
    
        print('Finished running foldseek screening')

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    if SCREEN == "fatcat":
        files = subprocess.check_output("ls", cwd=f"{_dir}/content/input", shell=True).decode().split("\n")[:-1]

        database_file = f"{_dir}/content/database/list_{DATABASE}.tab"
        with open(database_file, "r") as db:
            database_lines = [line.strip() for line in db.readlines()]  # Lê as linhas do arquivo de banco de dados

        fatcat_commands_file = f"{_dir}/fatcat_commands.txt"  # Arquivo para os comandos FATCATSearch.pl

        with open(fatcat_commands_file, "w") as file:
            for f in files:
                for line in database_lines:
                    subject = os.path.join(_dir, "content", "database", DATABASE, line + ".pdb")  # Constrói o caminho completo do subject
                    fatcat_command = [
                        f'{_dir}/content/programs/FATCAT-dist/FATCATMain/FATCAT',
                        '-p1',
                        f,
                        '-p2',
                        f"{line}.pdb",
                        "-b",
                        "-i1",
                        "./",
                        "-i2",
                        f"{_dir}/content/database/{DATABASE}"
                    ]
                    file.write(" ".join(fatcat_command) + "\n")  # Escreve o comando no arquivo

        # Executa os comandos usando o parallel
        with open(fatcat_commands_file, "r") as file:
            fatcat_output = subprocess.run(["parallel", "-j", str(num_processors)], stdin=file, capture_output=True, text=True).stdout

        # Processar os resultados individualmente para cada estrutura de entrada
        for f in files:
            tab_file = f"{_dir}/content/result/screening/{f}.tab.fmt"

            sort_command = ["sort", "-k11nr"]
            head_command = ["head", "-n", str(HEADN)]

            format_command = [
                "perl",
                f"{_dir}/content/programs/remolog/scripts/format_result_FATCAT.pl",
                "-",
                f"{_dir}/content/programs/remolog/data/maxScore_fatcat.tab"
            ]
            
            # Filtrar apenas os resultados para a estrutura atual
            grep_output = subprocess.check_output(["grep", f], input=fatcat_output.encode()).decode()

            # Ordenar e pegar os melhores resultados
            sorted_output = subprocess.check_output(sort_command, input=grep_output.encode()).decode()
            head_output = subprocess.check_output(head_command, input=sorted_output.encode()).decode()

            # Formatar a saída e salvar no arquivo tab_file
            formatted_output = subprocess.check_output(format_command, input=head_output.encode()).decode()
            with open(tab_file, "w") as file:
                file.write(formatted_output)

        fmt_files = glob.glob(f"{_dir}/content/result/screening/*.fmt")
        with open(f"{_dir}/content/result/fatcat_formatted.tab", "w") as output_file:
            for fmt_file in fmt_files:
                with open(fmt_file, "r") as input_file:
                    output_file.write(input_file.read())

        print('Finished running fatcat screening')

        end_time = time.time()  # Registrar o tempo de término da execução
        execution_time = end_time - start_time  # Calcular o tempo total de execução
        print(f"Total execution time: {execution_time} seconds")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    if SCREEN == "tmalign":
        start_time = time.time()
        files = subprocess.check_output("ls", cwd=f"{_dir}/content/input", shell=True).decode().split("\n")[:-1]

        tmalign_commands_file = f"{_dir}/tmalign_commands.txt"  # Arquivo para os comandos do TMalign

        # Ler a lista de nomes de proteínas do banco de dados
        database_file = f"{_dir}/content/database/list_{DATABASE}.tab"
        with open(database_file, "r") as db:
            database_lines = [line.strip() for line in db.readlines()]  # Lê as linhas do arquivo de banco de dados

        # Adicionar comandos TM-align e parser_cmd para cada proteína do banco de dados no arquivo externo
        with open(tmalign_commands_file, "w") as file:
            for f in files:
                for subject in database_lines:
                    tmalign_input_file = f"{_dir}/content/input/{f}"
                    tmalign_subject_file = f"{_dir}/content/database/{DATABASE}/{subject}.pdb"
                    parser_output_file = f"{_dir}/content/result/screening/{f}.tab"  # Arquivo para salvar a saída do parser

                    # Concatenar os comandos com o operador de pipeline "|"
                    cmd = f"{_dir}/content/bin/TMalign {tmalign_input_file} {tmalign_subject_file} | perl {_dir}/content/programs/remolog/scripts/parser_TMalign.pl - >> {parser_output_file}\n"
                    file.write(cmd)

        # Executa os comandos usando o parallel
        with open(tmalign_commands_file, "r") as file:
            parser_output = subprocess.check_output(["parallel", "-j", str(num_processors)], stdin=file).decode()
        tab_file = f"{_dir}/content/result/screening/" + f + ".tab"
        with open(tab_file, "a") as file:
            file.write(parser_output)
        
        for f in files:
            parser_output_file = f"{_dir}/content/result/screening/{f}.tab"
            sort_output = subprocess.check_output(["sort", "-k3nr", str(parser_output_file)]).decode()
            grep_output = subprocess.check_output(["grep", "-e", f], input=sort_output.encode()).decode()
            head_output = subprocess.check_output(["head", "-n", str(HEADN)], input=grep_output.encode()).decode()

            # output_file_path = f"{parser_output_file}.fmt"
            # with open(output_file_path, "w") as fmt_file:
            #     lines = head_output.split("\n")
            #     for line in lines:
            #         columns = line.split()
            #         if len(columns) > 1:
            #             columns[1] = os.path.splitext(columns[1])[0]  # Remover a extensão .pdb
            #         fmt_file.write("\t".join(columns) + "\n")

            
            output_file_path = f"{parser_output_file}.fmt"
            with open(output_file_path, "w") as fmt_file:
                lines = head_output.strip().split("\n")  # Remover espaços em branco no início e no final e, em seguida, dividir por linhas
                for line in lines:
                    if line.strip():  # Verificar se a linha não está vazia após remover espaços em branco
                        columns = line.split()
                        if len(columns) > 1:
                            columns[1] = os.path.splitext(columns[1])[0]  # Remover a extensão .pdb
                        fmt_file.write("\t".join(columns) + "\n")
        
        # Juntar os arquivos formatados em um único arquivo final .tab
        fmt_files = [file for file in os.listdir(f"{_dir}/content/result/screening") if file.endswith(".fmt")]
        with open(f"{_dir}/content/result/tmalign_formatted.tab", "w") as output_file:
            for fmt_file in fmt_files:
                with open(os.path.join(f"{_dir}/content/result/screening", fmt_file), "r") as input_file:
                    lines = input_file.readlines()
                    output_file.write("".join(lines))

        if os.path.isfile(tmalign_commands_file):
            os.remove(tmalign_commands_file)
       
        end_time = time.time()  # Registrar o tempo de término da execução
        execution_time = end_time - start_time  # Calcular o tempo total de execução
        print(f"Total execution time of tmalign: {execution_time} seconds")

        print('Finished running tmalign screening')

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    if SCREEN != "fatcat":
        start_time_fatcat = time.time()
        fatcatFile = f'{_dir}/content/result/fatcat_formatted.tab'
    
        if os.path.isfile(fatcatFile):
            os.remove(fatcatFile)
        fatcat_commands_file = f"{_dir}/fatcat_commands.txt"
        # Cria um arquivo temporário para armazenar os comandos do FatCat
        with open(fatcat_commands_file, 'w') as file:
            for f in os.listdir('.'):
                with open(f'{_dir}/content/result/screening/{f}.tab.fmt') as fatcatFileInputFile:
                    for l in fatcatFileInputFile.readlines():
                        columns = l.split('\t')
                        if len(columns) > 1:
                            l = columns[1].strip()
                            l = l + '.pdb'
                     
                        
                        cmd = f"{_dir}/content/programs/FATCAT-dist/FATCATMain/FATCAT -p1 {f} -il ./ -p2 {l} -i2 {_dir}/content/database/{DATABASE} -b | perl {_dir}/content/programs/remolog/scripts/format_result_FATCAT.pl - {_dir}/content/programs/remolog/data/maxScore_fatcat.tab >> {fatcatFile}\n"
                        file.write(cmd)

                        
    
        # Executa os comandos do FatCat em paralelo usando o arquivo temporário
        with open(fatcat_commands_file, 'r') as file:
            subprocess.run(["parallel", "-j", str(num_processors)], stdin=file, capture_output=True, text=True).stdout
    
      
        if os.path.isfile(f'{_dir}/fatcat_commands.txt'):
            os.remove(f'{_dir}/fatcat_commands.txt')
    
        if os.path.isfile(f'{_dir}/content/result/tmpfatcatfile'):
            os.remove(f'{_dir}/content/result/tmpfatcatfile')
    
        
        end_time_fatcat = time.time()  # Registrar o tempo de término da execução
        execution_time_fatcat = end_time_fatcat - start_time_fatcat  # Calcular o tempo total de execução
        print('Finished running fatcat')
        print(f"Total execution time of !fatcat: {execution_time_fatcat} seconds")

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    if SCREEN != "tmalign":
        start_tmalign_time = time.time()
        tmalignFile = f"{_dir}/content/result/tmalign_formatted.tab"
        tmalignOutput = f"{_dir}/content/result/tmalign.tab"
        tmalign_commands_file = f"{_dir}/tmalign_commands.txt"
        if os.path.isfile(tmalignFile):
            os.remove(tmalignFile)

        # Criar uma lista para armazenar os comandos TMalign
        tmalign_commands = []

        # Loop pelos arquivos .tab.fmt para gerar os comandos TMalign e armazená-los na lista
        for f in os.listdir('.'):
            with open(f'{_dir}/content/result/screening/{f}.tab.fmt') as tmalignInputFile:
                for line in tmalignInputFile.readlines():
                    columns = line.split('\t')
                    if len(columns) > 1:
                        db_protein = columns[1].strip() + '.pdb'
                        tabFile = f"{f}"
                        # Adicionar o comando TMalign à lista
                        tmalign_commands.append(f"{_dir}/content/bin/TMalign {tabFile} {_dir}/content/database/{DATABASE}/{db_protein} | perl {_dir}/content/programs/remolog/scripts/parser_TMalign.pl - >> {tmalignOutput}\n")

        # Salvar a lista de comandos em um arquivo
        with open(tmalign_commands_file, 'w') as file:
            file.write("\n".join(tmalign_commands))

        # Executar os comandos usando o parallel
        subprocess.run(["parallel", "-j", str(num_processors)], stdin=open(tmalign_commands_file), text=True)

        # Ler a saída do parallel sem aplicar formatação para remover a extensão .pdb da segunda coluna
        with open(tmalignOutput, "r") as parallel_output_file:
            parser_output = parallel_output_file.read()

        # Salvar o resultado formatado no arquivo tmalign_formatted.tab, removendo a extensão .pdb apenas da segunda coluna
        formatted_output_lines = []
        for line in parser_output.splitlines():
            columns = line.split('\t')
            if len(columns) > 1:
                columns[1] = columns[1].replace(".pdb", "")
                formatted_output_lines.append("\t".join(columns))

        formatted_output = "\n".join(formatted_output_lines)

        with open(tmalignFile, "w") as output_file:
            output_file.write(formatted_output)

        # Remover o arquivo de comandos do parallel
        if os.path.isfile(tmalign_commands_file):
            os.remove(tmalign_commands_file)

        if os.path.isfile(tmalignOutput):
            os.remove(tmalignOutput)

        end_tmalign_time = time.time()  # Registrar o tempo de término da execução
        execution_tmalign_time = end_tmalign_time - start_tmalign_time  # Calcular o tempo total de execução
        print('Finished running !tmalign')
        print(f"Total execution time of !tmalign: {execution_tmalign_time} seconds")


#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    start_lovoalign_time = time.time()

    if SCREEN != "lovoalign":
        lovoalignFile = f'{_dir}/content/result/lovoalign_formatted.tab'
        tempLovoalignFile = f'{_dir}/content/result/tempLovoalignFile.tab'
    
        if os.path.isfile(lovoalignFile):
            os.remove(lovoalignFile)
    
        lovoalign_commands_file = f"{_dir}/lovoalign_commands.txt"
        # Criar arquivo externo para armazenar os comandos do lovoalign
        with open(lovoalign_commands_file, 'w') as file:
            for f in os.listdir('.'):
                with open(f'{_dir}/content/result/screening/{f}.tab.fmt') as lovoalignInputFile:
                    for l in lovoalignInputFile.readlines():
                        columns = l.split('\t')
                        if len(columns) > 1:
                            l = columns[1].strip()
    
                                          
                        parser_output_file = f"{_dir}/content/result/lovotemp.tab"
                        cmd = f"{_dir}/content/bin/lovoalign -p1 {_dir}/content/input/{f} -p2 {_dir}/content/database/{DATABASE}/{l}.pdb | perl {_dir}/content/programs/remolog/scripts/parser_lovoalign.pl - >> {lovoalignFile}\n"
                        file.write(cmd)
    
        # Executar os comandos do lovoalign em paralelo usando o arquivo externo
        with open(lovoalign_commands_file, 'r') as file:
            subprocess.run(["parallel", "-j", str(num_processors)], stdin=file, capture_output=True, text=True).stdout
        
        # with open(parser_output_file, 'w') as lovoalignOutputFile:
        #     parserProcess = subprocess.Popen(parserCommand, stdout=subprocess.PIPE)
        #     parserOutput, _ = parserProcess.communicate()
        #     parserOutput = parserOutput.decode()
        #     lines = parserOutput.split("\n")
        #     for line in lines:
        #         columns = line.split("\t")
        #         if len(columns) > 1:
        #             columns[1] = os.path.splitext(columns[1])[0]
        #         formatted_line = "\t".join(columns)
        #         if formatted_line.strip():
        #             lovoalignOutputFile.write(formatted_line + "\n")
    
        if os.path.isfile(lovoalign_commands_file):
            os.remove(lovoalign_commands_file)
    
         
        end_lovoalign_time = time.time()  # Registrar o tempo de término da execução
        execution_lovoalign_time = end_lovoalign_time - start_lovoalign_time  # Calcular o tempo total de execução
        print(f"Total execution time of lovoalign: {execution_lovoalign_time} seconds")

        print('Finished running lovoalign')

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_

    result_dir = f"{_dir}/content/result"
    fatcat_file = f"{result_dir}/fatcat_formatted.tab"
    tmalign_file = f"{result_dir}/tmalign_formatted.tab"
    lovoalign_file = f"{result_dir}/lovoalign_formatted.tab"
    result_file = f"{result_dir}/result.tab"
   
    command1 = [
        "perl",
        f"{_dir}/content/programs/remolog/scripts/join_table.pl",
        fatcat_file,
        tmalign_file,
    ]
    processtest = subprocess.Popen(command1, stdout=subprocess.PIPE)
    tempCommand1File, _ = processtest.communicate()
                    
    with open(f'{_dir}/content/result/tempCommand1File', 'a') as fatcatOutputFile:
        fatcatOutputFile.write(tempCommand1File.decode())
        

    command2 = [
        "perl",
        f"{_dir}/content/programs/remolog/scripts/join_table.pl",
        f'{_dir}/content/result/tempCommand1File',
        lovoalign_file,
    ]

    processtest2 = subprocess.Popen(command2, stdout=subprocess.PIPE)
    tempCommand2File, _ = processtest2.communicate()
                
    with open(f'{_dir}/content/result/tempCommand2File', 'a') as fatcatOutputFile2:
        fatcatOutputFile2.write(tempCommand2File.decode())
    
    # Executar command3 usando o arquivo temporário como entrada
    command3 = [
        "perl",
        f"{_dir}/content/programs/remolog/scripts/add_scope_class.pl",
        f'{_dir}/content/result/tempCommand2File',
        ANNOT,
    ]

    # Redirecionar a saída para o arquivo result.tab
    with open(result_file, "w") as output_file:
        process3 = subprocess.Popen(
            command3, stdin=subprocess.PIPE, stdout=output_file, stderr=subprocess.PIPE
            )
        process3.wait()

    # # Remover arquivos temporarios
    # if os.path.isfile(f'{_dir}/content/result/tempCommand1File'):
    #     os.remove(f'{_dir}/content/result/tempCommand1File')
    # if os.path.isfile(f'{_dir}/content/result/tempCommand2File'):
    #     os.remove(f'{_dir}/content/result/tempCommand2File')
    # # Remover arquivos intermediários
    # if os.path.isfile(fatcatFile):
    #     os.remove(fatcatFile)
    # if os.path.isfile(tmalignFile):
    #     os.remove(tmalignFile)
    # if os.path.isfile(lovoalignFile):
    #     os.remove(lovoalignFile)
    
    print('Your task has been completed!')
    
    end_alltask_time = time.time()  # Registrar o tempo de término da execução
    execution_alltask_time = end_alltask_time - start_alltask_time  # Calcular o tempo total de execução
    print(f"Total execution time of your task: {execution_alltask_time} seconds")
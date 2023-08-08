import os
import subprocess

def create_database(_dir, database):
    print('*'*64, 'inside create_database')
    print(os.getcwd())
    os.chdir(_dir+'/content')

    DATABASE = database  # define o nome do banco de dados

    # cria o diretório para o banco de dados, se não existir
    # if not os.path.exists('./database'):
    #     os.mkdir('database')
    if not os.path.exists('database'):
        os.system('mkdir database')

    # os.chdir('database')
    os.chdir(_dir+'/content/database')

    # if not os.path.exists(DATABASE):
        # os.mkdir(DATABASE)
    if not os.path.exists(os.environ['DATABASE']):
        os.system(f'mkdir {DATABASE}')

    # os.chdir(DATABASE)
    os.chdir(os.environ['DATABASE'])
    
    # if not os.path.exists('../list_' + DATABASE + '.tab'): #inserir not
    #     if DATABASE == 'scope40':
    #         subprocess.run(['wget', 'https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-40-2.08.tgz'])
    #         subprocess.run(['tar', '-zxf', 'pdbstyle-sel-gs-bib-40-2.08.tgz'])
    #         subprocess.run(['find', 'pdbstyle-2.08', '-type', 'f', '-name', '*.ent', '-exec', 'mv', '{}', './', ';'])
    #         subprocess.run(['find . -name "*pdbstyle*" -exec rm -rf {} \;', '-execdir rm -rf {} \;'], shell=True)
    if not os.path.exists(f'{_dir}/content/database/list_{DATABASE}.tab'):
        if os.environ['DATABASE'] == 'scope40':
            os.system('wget https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-40-2.08.tgz')
            os.system('tar -zxf pdbstyle-sel-gs-bib-40-2.08.tgz')
            os.system('mv pdbstyle-2.08/*/*.ent .')
            os.system('rm -rf pdbstyle*')
            # subprocess.run(['find', 'pdbstyle-2.08', '-type', 'f', '-name', '*.ent', '-exec', 'mv', '{}', './', ';'])
            # subprocess.run(['find . -name "*pdbstyle*" -exec rm -rf {} \;', '-execdir rm -rf {} \;'], shell=True)

        # elif DATABASE == 'scope95':
        #     subprocess.run(['wget', 'https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-95-2.08.tgz'])
        #     subprocess.run(['tar', '-zxf', 'pdbstyle-sel-gs-bib-95-2.08.tgz'])
        #     subprocess.run(['find', 'pdbstyle-2.08', '-type', 'f', '-name', '*.ent', '-exec', 'mv', '{}', './', ';'])
        #     subprocess.run(['find . -name "*pdbstyle*" -exec rm -rf {} \;', '-execdir rm -rf {} \;'], shell=True)
        elif os.environ['DATABASE'] == 'scope95':
            os.system('wget https://scop.berkeley.edu/downloads/pdbstyle/pdbstyle-sel-gs-bib-95-2.08.tgz')
            os.system('tar -zxf pdbstyle-sel-gs-bib-95-2.08.tgz')
            os.system('mv pdbstyle-2.08/*/*.ent .')
            os.system('rm -rf pdbstyle*')
    
        # renomeia os arquivos .ent para .pdb
        for file in os.listdir():
            if file.endswith('.ent'):
                os.rename(file, file[:-4] + '.pdb')
        # os.chdir(_dir+'../list_$DATABASE.tab')
        # os.system('ls *.ent > ../list_$DATABASE.tab')
            
        with open('../list_' + DATABASE + '.tab', 'w') as f:
            for file in os.listdir():
                if file.endswith('.pdb'):
                    file_name = os.path.splitext(file)[0]  # Obtém o nome do arquivo sem a extensão
                    f.write(file_name + '\n')
        # os.system('for i in *.ent; do mv $i $i.pdb; done')
           
    # creating foldseek database
    os.chdir(_dir+'/content/foldseek_data')

    if not os.path.exists(f"fs_{database}"):
        os.system(f"{_dir}/content/bin/foldseek createdb {_dir}/content/database/{database} fs_{database}")
        
    os.chdir(_dir+'/content')
    os.system('if [ -d result ]; then rm -rf result ; fi')
    os.system('mkdir result;')
    os.system('mkdir result/screening;')
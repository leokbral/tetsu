import os
# import keyboard
#import subprocess

def get_dependencies(_dir):
    print('*'*64, 'inside get_dependencies')

    # Download e instalação do TMalign
    os.chdir(_dir+'/content/bin')
    if not os.path.exists('TMalign'):
        print('&'*64, 'Download e instalação do TMalign')
        os.system('wget "https://zhanggroup.org/TM-align/TMalign.cpp"')
        os.system('g++ -static -O3 -ffast-math -lm -o TMalign TMalign.cpp')
        
    # download and install FATCAT
    os.chdir(_dir+'/content/programs')
    if not os.path.exists('FATCAT-dist'):
        print('&'*64, 'download and install FATCAT')
        #subprocess.run(["git", "clone", "https://github.com/GodzikLab/FATCAT-dist.git"])
        #subprocess.run(['bash', '-c', 'cd FATCAT-dist; ./Install']) 
        os.system("if [ ! -d FATCAT-dist ]; then \
             git clone https://github.com/leokbral/FATCAT-dist.git; \
             cd FATCAT-dist/; ./Install; \
           fi")
        
    # download e instalação do lovoalign
    os.chdir(_dir+'/content/bin')
    if not os.path.exists('lovoalign-22.0.0'):
        #subprocess.run(['wget', 'https://github.com/m3g/lovoalign/archive/refs/tags/22.0.0.tar.gz'])
        #subprocess.run(['tar', '-xzf', '22.0.0.tar.gz'])
        #os.chdir('lovoalign-22.0.0/src')
        #subprocess.run(['make'])
        #subprocess.run(['cp', '../bin/lovoalign', './content'])
        #subprocess.run(['bash', '-c', 'cd ../../../programs && if [ ! -d remolog ]; then git clone https://github.com/tetsufmbio/remolog.git; fi'])
        #os.chdir('../../')
        print('&'*64, 'download e instalação do lovoalign')
        os.system('wget "https://github.com/m3g/lovoalign/archive/refs/tags/22.0.0.tar.gz"')
        os.system('tar -xzf 22.0.0.tar.gz')
        os.system('cd lovoalign-22.0.0/src && make')
        os.system('cp lovoalign-22.0.0/bin/lovoalign ./')
        
        
    # Download some scripts and model
    os.chdir(_dir+'/content/programs')
    if not os.path.exists('remolog'):
        #subprocess.run(["git", "clone", "https://github.com/tetsufmbio/remolog.git"])
        print('&'*64, 'Download some scripts and model')
        os.system('git clone https://github.com/leokbral/remolog.git')
    
    os.chdir(_dir+'/content/programs')
    # if not os.path.exists('foldseek'):
        #subprocess.run(['wget', 'https://mmseqs.com/foldseek/foldseek-linux-sse2.tar.gz'])
        #subprocess.run(['tar', '-xzf', 'foldseek-linux-sse2.tar.gz'])
        #subprocess.run(['cp', 'foldseek/bin/foldseek', '../bin/'])
        # print('&'*64, 'Download foldseek (o retorno)')
        # os.system('wget https://mmseqs.com/foldseek/foldseek-linux-sse2.tar.gz')
        # os.system('tar xvzf foldseek-linux-sse2.tar.gz')
        # os.system('cp foldseek/bin/foldseek ../bin')
    
    
    # if not os.path.exists('foldseek'):
    while not os.path.exists('foldseek'):
        # while True:
        try:
            print('&'*64, 'Download foldseek (o retorno)')
            # print('Pressione q para sair do loop...')
            os.system('wget https://mmseqs.com/foldseek/foldseek-linux-sse2.tar.gz')
            os.system('tar xvzf foldseek-linux-sse2.tar.gz')
            os.system('cp foldseek/bin/foldseek ../bin')
            # break  # Sair do loop se o código for executado com sucesso
        except:# Exception as e:
            # if "Unable to establish SSL connection." in str(e):
            # print(f"Erro DOIDO: {str(e)}")
            print('Tentando novamente...')
                 # # Aguardar até que a tecla 'q' seja pressionada para sair do loop
            # if keyboard.is_pressed('q'):
            #     print("Saindo do loop...")
            #     break

    print("Processo concluído com sucesso!")


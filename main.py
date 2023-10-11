import click
import os
import subprocess
import shutil

from get_dependencies import get_dependencies
from create_database import create_database
from all_in_one import running_proteins_structure


# identifica a pasta atual do script
_dir = os.path.dirname(os.path.abspath(__file__))

# define a pasta atual como o diretório de trabalho padrão
os.chdir(_dir)
if not os.path.exists(f'{_dir}/content/tmpDir'):
    os.makedirs(f'{_dir}/content/tmpDir')

@click.command()
@click.option('--screening', prompt='Choose an option for Screening software',
              type=click.Choice(['foldseek', 'tmalign', 'fatcat']),
              default='foldseek',
              help='Software to be used to retrieve structurally similar proteins in SCOPe database.')

@click.option('--num_analyze', prompt='Choose an option for most similar protein structures to be analyzed',
              type=int, default=200,
              help='Number of most similar protein structures to analyze.')

@click.option('--database', prompt='Choose an option for protein structure database to be used',
              type=click.Choice(['scope40', 'scope95']), default='scope40',
              help='Protein structure database to be used.')

@click.option('--job_name', prompt='Choose the jobName',
              type=str, default='remolog_final_result',
              help='Prefix for the final output name.')

@click.option('--num_processors', prompt='Choose the number of processors for parallel processing',
              type=int, default=12,
              help='Number of processors to use for parallel processing.')

# @click.option('--temp-dir', prompt='Enter the temporary directory path (leave blank for no temporary directory)',
#               type=click.Path(file_okay=False, writable=True, resolve_path=True), default='',
#               help='Directory to store temporary files.')

@click.option('--temp_dir', type=click.Path(exists=True), required=False, default=f'{_dir}/content/tmpDir')
@click.argument('input_dir', type=click.Path(exists=True), required=False)
def main(screening, num_analyze, database, job_name, num_processors, temp_dir, input_dir):
    if not os.path.exists(temp_dir):
       os.makedirs(temp_dir) 

    if input_dir:
        tmp_dir_path = os.environ.get('TMP_DIR_PATH', '')
        process_files(screening, num_analyze, database, job_name, num_processors, temp_dir, os.environ.get('TMP_DIR_PATH', ''), input_dir, tmp_dir_path)
    else:
        prompt_inputs(screening, num_analyze, database, job_name, num_processors, temp_dir)


def process_files(screening, num_analyze, database, job_name, num_processors, temp_dir, temp_dir_input, input_dir, tmp_dir_path):
    print(f'The chosen option was {screening}')
    print(f"Number of protein structures to analyze: {num_analyze}")
    print(f'Protein structure database to be used: {database}')
    print(f'jobName to be used: {job_name}')
    print(f'Number of processors choosed: {num_processors}')
    print(f'Input directory: {input_dir}')

    os.environ['FATCAT'] = f'{_dir}/content/programs/FATCAT-dist'
    os.environ['PATH'] += f':{_dir}/content/programs/FATCAT-dist/FATCATMain'
    os.environ['HEADN'] = str(num_analyze)
    os.environ['SCREEN'] = screening
    os.environ['DATABASE'] = database
    annot = ''

    if database == "scope40":
        os.environ['ANNOT'] = f"{_dir}/content/programs/remolog/data/dir.cla.scope.2.08-stable_filtered40.txt"
        annot = os.environ['ANNOT']
    elif database == "scope95":
        os.environ['ANNOT'] = f"{_dir}/content/programs/remolog/data/dir.cla.scope.2.08-stable_filtered95.txt"
        annot = os.environ['ANNOT']

    if not os.path.exists(f'{_dir}/content'):
        os.makedirs(f'{_dir}/content')
        os.makedirs(f'{_dir}/content/input/')
        # os.makedirs(f'{_dir}/content/tmpDir')####

    input_dest = os.path.join(_dir, 'content', 'input/')####
    if input_dir:
        for file_name in os.listdir(input_dir):
            file_path = os.path.join(input_dir, file_name)
            if os.path.isfile(file_path):
                shutil.copy(file_path, input_dest)
    
    # if temp_dir:
    #     os.environ['TMP_DIR'] = temp_dir
    # else:
    #     os.environ['TMP_DIR'] = f'{_dir}/content/tmp/'

    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d bin ]; then mkdir bin; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d programs ]; then mkdir programs; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d view ]; then mkdir view; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d foldseek_data ]; then mkdir foldseek_data; fi'])

    get_dependencies(_dir)
    create_database(_dir, database)
    running_proteins_structure(_dir, screening, database, num_analyze, annot, num_processors, temp_dir_input)
   


def prompt_inputs(screening, num_analyze, database, job_name, num_processors, temp_dir):
    print(f'The chosen option was {screening}')
    print(f"Number of protein structures to analyze: {num_analyze}")
    print(f'Protein structure database to be used: {database}')
    print(f'jobName to be used: {job_name}')
    print(f'Number of processors choosed: {num_processors}')
    
    while True:
        temp_dir_input = click.prompt('Enter the temporary directory path (leave blank for no change)',
                                      default=temp_dir, type=click.Path(file_okay=False, writable=True, resolve_path=True))
        
        if not temp_dir_input:
            break
    
        if os.path.exists(temp_dir_input) and os.path.isdir(temp_dir_input):
            os.environ['TMP_DIR_PATH'] = temp_dir_input  # Set the environment variable
            break
        else:
            create_dir = click.confirm(f'The path "{temp_dir_input}" is not a valid directory. Would you like to create this directory?', default=False)
            if create_dir:
                try:
                    os.makedirs(temp_dir_input)
                    os.environ['TMP_DIR_PATH'] = temp_dir_input  # Set the environment variable
                    break
                except Exception as e:
                    print(f'Error creating the directory: {e}')
            else:
                print('Enter a valid path or leave it blank.')

    while True:
        input_dir = click.prompt('Enter the input directory path', type=click.Path(exists=True))
        if os.path.isdir(input_dir):
            break
        else:
            print(f'Error: Path "{input_dir}" is not a directory or does not exist.')

    process_files(screening, num_analyze, database, job_name, num_processors, temp_dir, temp_dir_input, input_dir, os.environ.get('TMP_DIR_PATH', ''))


if __name__ == '__main__':
    main()

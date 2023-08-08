import click
import os
import subprocess
import shutil

from get_dependencies import get_dependencies
from create_database import create_database
# from running_protein_structure_alignment import running_proteins_structure #COMENTAR
# from parallel_test import parallel_test
# from new import run_script
# from paral import run_script
# from rpsa import running_proteins_structure
# from fatcat_function import running_proteins_structure#, running_fatcat_parallel #DESCOMENTAR
# from tmalign_function import running_proteins_structure
# from lovoalign_function import running_proteins_structure
from all_in_one import running_proteins_structure
# from fatcat_function import *

# identifica a pasta atual do script
_dir = os.path.dirname(os.path.abspath(__file__))

# define a pasta atual como o diretório de trabalho padrão
os.chdir(_dir)

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

@click.argument('input_dir', type=click.Path(exists=True), required=False)
def main(screening, num_analyze, database, job_name, num_processors, input_dir):
    if input_dir:
        process_files(screening, num_analyze, database, job_name, num_processors, input_dir)
    else:
        prompt_inputs(screening, num_analyze, database, job_name, num_processors)


def process_files(screening, num_analyze, database, job_name, num_processors, input_dir):
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
        os.makedirs(f'{_dir}/content/input/')####

    input_dest = os.path.join(_dir, 'content', 'input/')####
    if input_dir:
        for file_name in os.listdir(input_dir):
            file_path = os.path.join(input_dir, file_name)
            if os.path.isfile(file_path):
                shutil.copy(file_path, input_dest)

    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d bin ]; then mkdir bin; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d programs ]; then mkdir programs; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d view ]; then mkdir view; fi'])
    subprocess.run(['bash', '-c', 'cd ./content/ && if [ ! -d foldseek_data ]; then mkdir foldseek_data; fi'])

    get_dependencies(_dir)
    create_database(_dir, database)
    running_proteins_structure(_dir, screening, database, num_analyze, annot, num_processors)
    # running_fatcat_parallel(_dir, screening, database, num_analyze, annot)
    # fatcattest(_dir, screening, database, num_analyze, annot)
    # parallel_test(_dir, screening, database, num_analyze, annot)
    # run_script(_dir, screening, database, num_analyze, annot)



def prompt_inputs(screening, num_analyze, database, job_name, num_processors):
    print(f'The chosen option was {screening}')
    print(f"Number of protein structures to analyze: {num_analyze}")
    print(f'Protein structure database to be used: {database}')
    print(f'jobName to be used: {job_name}')
    print(f'Number of processors choosed: {num_processors}')

    while True:
        input_dir = click.prompt('Enter the input directory path', type=click.Path(exists=True))
        if os.path.isdir(input_dir):
            break
        else:
            print(f'Error: Path "{input_dir}" is not a directory or does not exist.')

    process_files(screening, num_analyze, database, job_name, input_dir)


if __name__ == '__main__':
    main()

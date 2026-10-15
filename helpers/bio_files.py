from pathlib import Path

def prepare_fasta(input_fasta, output_fasta):

    input_fasta_path = Path(input_fasta)

    file = open(input_fasta_path, 'r')

    if output_fasta is not None:
        output_fasta_path = Path(output_fasta)
    else:
        input_dir, input_name = input_fasta_path.parent, input_fasta_path.stem
        output_fasta_path = Path(input_dir) / f'{input_name}_uniline.fasta'

    output_file = open(output_fasta_path, 'w')

    return file, output_file

def write_fasta(s, output_file):
    output_file.write(s+'\n')

def prepare_blast(input_file, output_file):

    input_file_path = Path(input_file)

    file = open(input_file_path, 'r')

    if output_file is not None:
            output_file_path = Path(output_file)
    else:
        input_dir, input_name = input_file_path.parent, input_file_path.stem
        output_file_path = Path(input_dir) / f'{input_name}_selected_blast_queries.txt'
    
        output_file = open(output_file_path, 'w')
    
    return file, output_file

def write_file(names, output_file):
    for name in names:
        output_file.write(name +'\n')
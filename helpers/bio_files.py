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
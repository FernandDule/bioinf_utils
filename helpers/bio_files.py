from pathlib import Path
from typing import TextIO

def prepare_fasta(input_fasta: str,
                  output_fasta: str | None = None) -> tuple[TextIO, TextIO]:
    '''
    Opens input fasta file for reading and prepares output file for writing.
    If output path is not provided, creates file in the input directory with '_uniline.fasta' suffix.

    Arguments:
        input_fasta: str - path to input fasta file.
        output_fasta (opt): str - path to output fasta file. Default = None.
    '''

    input_fasta_path = Path(input_fasta)

    file = open(input_fasta_path, 'r')

    if output_fasta is not None:
        output_fasta_path = Path(output_fasta)
    else:
        input_dir, input_name = input_fasta_path.parent, input_fasta_path.stem
        output_fasta_path = Path(input_dir) / f'{input_name}_uniline.fasta'

    output_file = open(output_fasta_path, 'w')

    return file, output_file

def write_fasta(s: str,
                output_file: TextIO) -> None:
    '''
    Writes a uniline fasta sequence string to the output file.

    Arguments:
        s: str - sequence string.
        output_file: TextIO - opened output file object.
    '''
    output_file.write(s+'\n')

def prepare_blast(input_file: str, 
                  output_file: str | None = None) -> tuple[TextIO, TextIO]:
    '''
    Opens input blast file for reading and prepares output file for writing.
    If output path is not provided, creates it in the input directory with '_selected_blast_queries.txt' suffix.

    Arguments:
        input_file: str - path to input blast results file.
        output_file (opt): str - path to output file. Default = None.
    '''

    input_file_path = Path(input_file)

    file = open(input_file_path, 'r')

    if output_file is not None:
        output_file_path = Path(output_file)
        output_file = open(output_file_path, 'w')
    else:
        input_dir, input_name = input_file_path.parent, input_file_path.stem
        output_file_path = Path(input_dir) / f'{input_name}_selected_blast_queries.txt'

        output_file = open(output_file_path, 'w')

    return file, output_file

def write_file(names: list[str],
               output_file: TextIO) -> None:
    '''
    Writes a list of protein names to the output file, each on a new line.

    Arguments:
        names: list[str] - list of protein names.
        output_file: TextIO - opened output file object.
    '''
    for name in names:
        output_file.write(name +'\n')
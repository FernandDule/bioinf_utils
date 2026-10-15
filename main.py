from pathlib import Path

import sys
sys.path.insert(0, "./helpers")

from fastq import open_fastq, gc_control, length_control, quality_control, write_fastq, parse_bounds, write_log, write_failed
from bio_files import prepare_fasta, write_fasta, prepare_blast, write_file


def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):

    log = open('./log.txt', 'w')

    output_file_path = Path(output_fastq)
    output_dir = Path('./filtered')
    output_dir.mkdir(exist_ok=True, parents=True)

    filtered = open(output_dir/output_file_path, 'w', encoding='ascii')

    min_gc, max_gc = parse_bounds(gc_bounds)
    min_length, max_length = parse_bounds(length_bounds)

    n_total = 0
    n_gc_fail = 0
    n_length_fail = 0
    n_quality_fail = 0

    file = open_fastq(input_fastq)

    while True:

        name = file.readline().rstrip('\n')
        if not name:
            break

        n_total += 1

        sequence = file.readline().rstrip('\n')
        file.readline()
        quality = file.readline().rstrip('\n')

        if not gc_control(sequence, min_gc, max_gc):
            n_gc_fail += 1
            write_failed(log, name, 'gc_control')
            continue

        if not length_control(sequence, min_length, max_length):
            n_length_fail += 1
            write_failed(log, name, 'length_control')
            continue

        if not quality_control(quality, quality_threshold):
            n_quality_fail += 1
            write_failed(log, name, 'quality_control')
            continue

        write_fastq(filtered, name, sequence, quality)

    write_log(log, input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold,
              n_total, n_gc_fail, n_length_fail, n_quality_fail)

    file.close()
    filtered.close()
    log.close()

def convert_multiline_fasta_to_oneline(input_fasta, output_fasta=None):

    file, output_file = prepare_fasta(input_fasta, output_fasta)

    count = 0
    s = ""

    while True:
        line = file.readline().rstrip('\n')

        if not line:
            break

        if '>' in line:
            if s:
                write_fasta(s, output_file)
                count += 1

            s = f'{line}\n'
        else:
            s += line

    if s:
        write_fasta(s, output_file)
        count += 1

    print(f'Total number of fasta-sequences: {count}')
    file.close()
    output_file.close()

def parse_blast_output(input_file, output_file=None):

    file, output_file = prepare_blast(input_file, output_file)

    n_queries = 0
    names = []

    while True:

            line = file.readline()

            if not line:
                break

            if 'Sequences producing significant alignments:' in line:
                n_queries += 1

                file.readline()
                file.readline()
                line = file.readline()
                accession = line.rstrip('\n').split()[-1]

            if 'Alignments:' in line:
                file.readline()
                name_line = file.readline()
                accession_line = file.readline()

                if '>' in name_line and (accession in accession_line):
                    name = name_line.lstrip('>').replace('MULTISPECIES: ', '').rstrip('\n')
                    names.append(name)


    write_file(sorted(names), output_file)
    print(f'Number of queries: {n_queries}')
    print(f'Proteins saved: {len(names)}')

    file.close()
    output_file.close()
    
from pathlib import Path

import sys
sys.path.insert(0, "./helpers")
from fastq import open_fastq, gc_control, length_control, quality_control, write_fastq

def filter_fastq(input_fastq, output_fastq, gc_bounds=(0, 100), length_bounds=(0, 2**32), quality_threshold=0):

    log = open('./log.txt', 'w')
    log.write(f'input_fastq file: {input_fastq} \n')
    log.write(f'output_fastq file: {output_fastq} \n')
    log.write(f'gc_bounds: {gc_bounds} \n')
    log.write(f'length_counds: {length_bounds} \n')
    log.write(f'quality_threshold (phred33): {quality_threshold} \n')

    output_file_path = Path(output_fastq)
    output_dir = output_file_path.parent
    output_dir.mkdir(exist_ok=True, parents=True)
    filtered = open(output_file_path, 'w', encoding='ascii')

    if len(gc_bounds) == 2:
        min_gc, max_gc = gc_bounds
    else:
        min_gc, max_gc = 0, gc_bounds

    if len(length_bounds) == 2:
        min_length, max_length = length_bounds
    else:
        min_length, max_length = 0, length_bounds

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
            log.write(f'{name}\n')
            log.write('Failed gc control \n')
            continue

        if not length_control(sequence, min_length, max_length):
            n_length_fail += 1
            log.write(f'{name}\n')
            log.write('Failed length control \n')
            continue

        if not quality_control(quality, quality_threshold):
            n_quality_fail += 1
            log.write(f'{name}\n')
            log.write('Failed quality control \n')
            continue

        write_fastq(filtered, name, sequence, quality)

    print(f'Initial number of sequences: {n_total} \n')
    print(f'Failed gc_control: {n_gc_fail} ({n_gc_fail/n_total*100:.2f}%) \n')
    print(f'Failed length_control: {n_length_fail} ({n_length_fail/n_total*100:.2f}%) \n')
    print(f'Failed quality_control: {n_quality_fail} ({n_quality_fail/n_total*100:.2f}%) \n')
    print(f'Total failed: {n_gc_fail + n_length_fail + n_quality_fail} ({(n_gc_fail + n_length_fail + n_quality_fail)/n_total:.2f}%)')
    print(f'Passed: {n_total - (n_gc_fail + n_length_fail + n_quality_fail)} ({(n_total - (n_gc_fail + n_length_fail + n_quality_fail))/n_total:.2f}%)')

    log.write(f'Initial number of sequences: {n_total} \n')
    log.write(f'Failed gc_control: {n_gc_fail} ({n_gc_fail/n_total*100:.2f}%) \n')
    log.write(f'Failed length_control: {n_length_fail} ({n_length_fail/n_total*100:.2f}%) \n')
    log.write(f'Failed quality_control: {n_quality_fail} ({n_quality_fail/n_total*100:.2f}%) \n')
    log.write(f'Total failed: {n_gc_fail + n_length_fail + n_quality_fail} ({(n_gc_fail + n_length_fail + n_quality_fail)/n_total:.2f}%)')
    log.write(f'Passed: {n_total - (n_gc_fail + n_length_fail + n_quality_fail)} ({(n_total - (n_gc_fail + n_length_fail + n_quality_fail))/n_total:.2f}%)')
    
    file.close()
    filtered.close()
    log.close()

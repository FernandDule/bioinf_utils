from pathlib import Path
from typing import TextIO

def parse_bounds(bounds: tuple[int | float, int | float] | int | float) -> list[int | float]:

    '''
    Parses bounds argument. If a single number is passed, sets min bound to 0.

    Arguments:
        bounds: tuple[int | float, int | float] | int | float - min and max bounds.
    '''

    if isinstance(bounds, tuple):
        bounds = [bounds[0], bounds[1]]
    else:
        bounds = [0, bounds]
    return bounds

def open_fastq(input_fastq: str) -> TextIO:

    '''
    Opens input fastq file for reading in ascii encoding.

    Arguments:
        input_fastq: str - path to input fastq file.
    '''

    input_fastq_path = Path(input_fastq)
    return open(input_fastq_path, encoding='ascii')

def gc_control(sequence: str,
               min_gc: int | float, max_gc: int | float) -> bool:

    '''
    Checks if GC content of the sequence is within specified bounds.

    Arguments:
        sequence: str - sequence.
        min_gc: int | float - minimum GC percentage.
        max_gc: int | float - maximum GC percentage.
    '''

    seq = sequence.upper()
    gc = (seq.count('C') + seq.count('G')) / len(seq)*100
    return min_gc <= gc <= max_gc

def length_control(sequence: str,
                   min_length: int,
                   max_length: int) -> bool:

    '''
    Checks if sequence length is within specified bounds.

    Arguments:
        sequence: str - sequence.
        min_length: int - minimum sequence length.
        max_length: int - maximum sequence length.
    '''
    return min_length <= len(sequence) <= max_length

def quality_control(quality: str,
                    quality_threshold: int | float) -> bool:
    '''
    Checks if mean phred33 quality score of the sequence is above threshold.

    Arguments:
        quality: str - quality string from fastq file.
        quality_threshold: int | float - minimum average quality score.
    '''

    sum_score = 0
    for ch in quality:
        score = ord(ch) - 33
        sum_score += score

    return sum_score/len(quality) >= quality_threshold

def write_fastq(filtered: TextIO,
                name: str,
                sequence: str,
                quality: str) -> None:
    '''
    Writes a fastq sequence to the output file.

    Arguments:
        filtered: TextIO - opened output file object.
        name: str - sequence header.
        sequence: str - sequence string.
        quality: str - quality string.
    '''

    filtered.write(f'{name}\n')
    filtered.write(f'{sequence}\n')
    filtered.write(f'+\n')
    filtered.write(f'{quality}\n')

def write_log(log: TextIO,
              input_fastq: str,
              output_fastq: str,
              gc_bounds: tuple[int, int] | int | float,
              length_bounds: tuple[int, int] | int | float,
              quality_threshold: int | float,
              n_total: int,
              n_gc_fail: int,
              n_length_fail: int,
              n_quality_fail: int) -> None:

    '''
    Writes filtering statistics to log file and prints them to console.

    Arguments:
        log: TextIO - opened log file object.
        input_fastq: str - path to input fastq file.
        output_fastq: str - path to output fastq file.
        gc_bounds: tuple[int, int] | int | float - GC content bounds used.
        length_bounds: tuple[int, int] | int | float - length bounds used.
        quality_threshold: int | float - quality threshold used.
        n_total: int - total number of sequences.
        n_gc_fail: int - number of sequences failed GC control.
        n_length_fail: int - number of sequences failed length control.
        n_quality_fail: int - number of sequences failed quality control.
    '''

    log.write(f'input_fastq file: {input_fastq} \n')
    log.write(f'output_fastq file: {output_fastq} \n')
    log.write(f'gc_bounds: {gc_bounds} \n')
    log.write(f'length_bounds: {length_bounds} \n')
    log.write(f'quality_threshold (phred33): {quality_threshold} \n')

    print(f'Initial number of sequences: {n_total} \n')
    print(f'Failed gc_control: {n_gc_fail} ({n_gc_fail/n_total*100:.2f}%) \n')
    print(f'Failed length_control: {n_length_fail} ({n_length_fail/n_total*100:.2f}%) \n')
    print(f'Failed quality_control: {n_quality_fail} ({n_quality_fail/n_total*100:.2f}%) \n')
    print(f'Total failed: {n_gc_fail + n_length_fail + n_quality_fail} ({(n_gc_fail + n_length_fail + n_quality_fail)/n_total*100:.2f}%)')
    print(f'Passed: {n_total - (n_gc_fail + n_length_fail + n_quality_fail)} ({(n_total - (n_gc_fail + n_length_fail + n_quality_fail))/n_total*100:.2f}%)')

    log.write(f'Initial number of sequences: {n_total} \n')
    log.write(f'Failed gc_control: {n_gc_fail} ({n_gc_fail/n_total*100:.2f}%) \n')
    log.write(f'Failed length_control: {n_length_fail} ({n_length_fail/n_total*100:.2f}%) \n')
    log.write(f'Failed quality_control: {n_quality_fail} ({n_quality_fail/n_total*100:.2f}%) \n')
    log.write(f'Total failed: {n_gc_fail + n_length_fail + n_quality_fail} ({(n_gc_fail + n_length_fail + n_quality_fail)/n_total*100:.2f}%)')
    log.write(f'Passed: {n_total - (n_gc_fail + n_length_fail + n_quality_fail)} ({(n_total - (n_gc_fail + n_length_fail + n_quality_fail))/n_total*100:.2f}%)')


def write_failed(log: TextIO,
                 name: str,
                 filter_name: str) -> None:
    '''
    Writes failed sequence header and filter name to log file.

    Arguments:
        log: TextIO - opened log file object.
        name: str - sequence header.
        filter_name: str - name of the failed filter.
    '''

    log.write(f'{name}\n')
    log.write(f'Failed {filter_name} \n')
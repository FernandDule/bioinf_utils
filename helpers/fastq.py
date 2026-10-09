from pathlib import Path

def parse_bounds(bounds):
    s = bounds.strip('()')
    parts = [p.strip() for p in s.split(',') if p.strip()]
    bounds = tuple(float(p) for p in parts) if len(parts) > 1 else (0, float(parts[0]))
    return bounds

def open_fastq(input_fastq):
    input_fastq_path = Path(input_fastq)
    return open(input_fastq_path, encoding='ascii')

def gc_control(sequence, min_gc, max_gc):

    seq = sequence.upper()
    gc = (seq.count('C') + seq.count('G')) / len(seq)*100
    return min_gc <= gc <= max_gc

def length_control(sequence, min_length, max_length):
    return min_length <= len(sequence) <= max_length

def quality_control(quality, quality_threshold):
   sum_score = 0
   for ch in quality:
    score = ord(ch) - 33
    sum_score += score

   return sum_score/len(quality) >= quality_threshold

def write_fastq(filtered, name, sequence, quality):
    filtered.write(f'{name}\n')
    filtered.write(f'{sequence}\n')
    filtered.write(f'+\n')
    filtered.write(f'{quality}\n')

def write_log(log, input_fastq, output_fastq, gc_bounds, length_bounds, quality_threshold, 
              n_total, n_gc_fail, n_length_fail, n_quality_fail):
    log.write(f'input_fastq file: {input_fastq} \n')
    log.write(f'output_fastq file: {output_fastq} \n')
    log.write(f'gc_bounds: {gc_bounds} \n')
    log.write(f'length_counds: {length_bounds} \n')
    log.write(f'quality_threshold (phred33): {quality_threshold} \n')

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


def write_failed(log, name, filter):
    log.write(f'{name}\n')
    log.write(f'Failed {filter} \n')
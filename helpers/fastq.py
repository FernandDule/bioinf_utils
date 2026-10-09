from pathlib import Path

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
   
   return sum_score/len(quality) > quality_threshold

def write_fastq(filtered, name, sequence, quality):
    filtered.write(f'{name}\n')
    filtered.write(f'{sequence}\n')
    filtered.write(f'+\n')
    filtered.write(f'{quality}\n')
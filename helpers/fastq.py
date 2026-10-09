def gc_control(sequence, min_gc, max_gc):
    
    seq = sequence.upper()
    gc = (seq.count('C') + seq.count('G')) / len(seq)
    return min_gc <= gc <= max_gc

def length_control(sequence, min_length, max_length):
    return min_length <= len(sequence) <= max_length

def quality_control(quality, quality_threshold):
   sum_score = 0
   for ch in quality:
    score = ord(ch) - 33
    sum_score += score
   
   return sum_score/len(quality) > quality_threshold

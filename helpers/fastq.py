def gc_control(sequence, min_gc, max_gc):
    
    seq = sequence.upper()
    gc = (seq.count('C') + seq.count('G')) / len(seq)
    return min_gc <= gc <= max_gc
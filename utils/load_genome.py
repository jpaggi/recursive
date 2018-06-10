def load_genome(genome):
    chrom = ''
    seq = []
    out = {}
    for line in genome:
        line = line.strip()
        if not line: continue
        if line[0] == '>':
            out[chrom] = ''.join(seq)
            seq = []
            chrom = line[1:].split(' ')[0]
        else:
            seq += [line]
    out[chrom] = ''.join(seq)
    return out

def complement(char):
    if char == 'A': return 'T'
    if char == 'C': return 'G'
    if char == 'G': return 'C'
    if char == 'T': return 'A'
    return char

def revcomp(seq):
    return ''.join(map(complement, seq[::-1]))
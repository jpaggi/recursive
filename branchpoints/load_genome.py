def load_genome(genome):
    chrom = ''
    seq = []
    out = {}
    for line in genome:
        line = line.strip()
        if line[0] == '>':
            out[chrom] = ''.join(seq)
            seq = []
            chrom = line[1:].split(' ')[0]
        else:
            seq += [line]
    out[chrom] = ''.join(seq)
    return out

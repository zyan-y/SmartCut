

def rev_comp(s):
    bp = {'A':'T', 'T':'A', 'G':'C', 'C':'G'}
    sc = ''.join(bp[i] for i in reversed(s))
    return sc


def idx_to_oligo(sequence, preds):
    oligos = []
    primers = []
    for i in range(len(preds)):
        if i == 0:
            oligo = sequence[:preds[i][1]]
        else:
            oligo = sequence[preds[i-1][0]:preds[i][1]]
        if (i%2) != 0:
            oligo = rev_comp(oligo)
        oligos.append(oligo)
        primers.append(sequence[preds[i][0]:preds[i][1]])
    oligos.append(rev_comp(sequence[preds[i][0]:]))
    return oligos, primers
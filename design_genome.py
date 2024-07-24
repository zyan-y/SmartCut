import os
import time
import pandas as pd
from Bio import SeqIO

from SmartCut.design_modules import assembly_design


def design_block(name, seq, savepath, expect_length=80, expect_overlap=18, flexibility=9, Tm_range=7):
    oligos, overlaps = assembly_design(seq, Tm_range, expect_length, expect_overlap, flexibility)
    overlaps.append('')
    results = [['idx', 'oligo', 'length of oligo', 'overlap', 'length of overlap']]
    for i in range(len(oligos)):
        results.append(['{}-{}'.format(name, i+1), oligos[i], len(oligos[i]), overlaps[i], len(overlaps[i])])
    results[-1][-1] = ''
    if savepath:
        results = pd.DataFrame(results)
        results.to_excel(savepath + name +'-oligos.xlsx', header=False, index=False)


def design_genome(genome, seq, block_length, block_overlap=0):
    seq = seq.upper()
    if len(seq) < block_length:
        print('Sequence length is too short!')
    else:
        num = round(len(seq)/block_length)
        blocks = []
        for n in range(num):
            BuildingBlock = seq[n*block_length: (n+1)*block_length+block_overlap]
            blocks.append([f'{genome}-block-{n+1}', BuildingBlock])
        if len(seq) > num*block_length+block_overlap:
            blocks[-1] = [f'{genome}-block-{num}', seq[(num-1)*block_length:]]
            
        print(f'Will design {len(blocks)} blocks for {genome}.')
        
        savepath = f'./genome_analysis/{genome}/'
        if not os.path.exists(savepath):
            os.mkdir(f'./genome_analysis/{genome}/')
        else:
            print('Path already exists! New oligos will over-write the folder.')
            
        for name, seq in blocks:
            design_block(name, seq, savepath)


if __name__ == "__main__":
    t1 = time.perf_counter()

    genome = 'Yeast_synV'

    seq = str(SeqIO.read(f'./genome_analysis/data/{genome}.fa', 'fasta').seq)
    block_length = 1000
    design_genome(genome, seq, block_length)
    
    t2 = time.perf_counter()
    print(f'Design {genome} cost {t2-t1}s')

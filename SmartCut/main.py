
import argparse
import time
import pandas as pd
import design_modules



def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('input_file', type=str, help='Path to the input sequence file')
    parser.add_argument('expect_length', type=int, help='Expected length')
    parser.add_argument('expect_overlap', type=int, help='Expected overlap')
    parser.add_argument('flexibility', type=int, help='Flexibility')
    parser.add_argument('Tm_range', type=int, help='Tm range')
    parser.add_argument('savepath', type=str, help='Path to save the results')
    parser.add_argument('diff', type=bool, help='Whether to predict synthetic difficulty')
    
    args = parser.parse_args()

    
    t1 = time.perf_counter()

    seqs = pd.read_excel(args.input_file, header=0, index_col=None, dtype='str').values[:,0:2]
    for name, seq in seqs:
        seq = seq.upper()
        if len(seq) > 3000:
            print('May be too long for PCA, still continue')
        oligos, overlaps = design_modules.assembly_design(seq, args.Tm_range, args.expect_length, args.expect_overlap, args.flexibility)

        results = [['idx', 'oligo', 'length of oligo', 'overlap', 'length of overlap']]
        for i in range(len(oligos)):
            results.append(['{}-{}'.format(name, i+1), oligos[i], len(oligos[i]), overlaps[i], len(overlaps[i])])
        results[-1][-1] = ''
        results = pd.DataFrame(results)
        results.to_excel(args.savepath + name+'.xlsx', header=None, index=None)

    t2 = time.perf_counter()    
    print('finish {}, use {:.3f}s'.format(args.input_file, t2-t1))


if __name__ == "__main__":
    main()
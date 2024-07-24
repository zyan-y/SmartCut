import sys
sys.path.append('./SmartCut')
import traceback
import numpy as np
import primer3
from model import predict_data
from utlis import *

step_value = 1
adjust_limit = 10


def cutoff(sequence, stride, expect_overlap):
    infor = [] # seq, start
    sequence = sequence[:-expect_overlap]

    num = round(len(sequence)/stride)
    oligo_num = num if num % 2 == 0 else num + 1
    oligo_length = int(len(sequence)/oligo_num)
    
    # adjust the length of the last
    excess = len(sequence) - oligo_length*oligo_num
    if excess <= 5:
        for i in range(1, oligo_num):
            infor.append([sequence[oligo_length*i:oligo_length*(i+1)], oligo_length*i])
    else:
        for i in range(1, excess):
            length = oligo_length + 1
            infor.append([sequence[length*i:length*(i+1)], length*i])
        n = length*excess
        for i in range(0, oligo_num - excess):
            infor.append([sequence[n + oligo_length*i: n + oligo_length*(i+1)], n + oligo_length*i])
    return infor


def data_processing(infor, domain, expect_overlap):
    data = [] 
    for i in range(len(infor)):
        seq, seq_start = infor[i][0], infor[i][1]
        candidates = []
        for j in range(domain//step_value):
            start = j * step_value
            input_seq = seq[start: start + expect_overlap]            
            candidates.append([input_seq, seq_start + start, seq_start + start + expect_overlap])
        data.append(candidates)
    return data


def predict(data):
    input_seqs = [j[0] for i in data for j in i] 
    scores = predict_data(input_seqs)
    scores = scores.reshape((len(data), -1))
    pred = []
    for i in range(len(data)):
        rank = np.argmax(scores[i])
        pred.append(data[i][rank][1:])
    return pred


def adjust_Tm(sequence, preds, Tm_range):

    Tms = []

    for pos in preds:
        primer = sequence[pos[0]: pos[1]]
        Tm = primer3.calc_tm(primer)
        Tms.append(Tm)

    retry_count = 0
    while np.ptp(Tms) > Tm_range*2:
        retry_count = retry_count + 1
        if retry_count > adjust_limit:
            break

        median = np.median(Tms)
        n,dist = 0,0
        for i in range(len(Tms)):
            distance = abs(Tms[i] - median)
            if distance > dist:
                n = i
                dist = distance
                
        limit = 3
        param = -1 if Tms[n] - median > 0 else 1
        if sequence[preds[n][1]-1] in 'CG':
            for i in range(1, limit+1):
                primer = sequence[preds[n][0]-i*param : preds[n][1]]
                Tm = primer3.calc_tm(primer)
                if abs(Tm - median) <= Tm_range:
                    break
            Tms[n] = Tm
            preds[n] = [preds[n][0]-i*param, preds[n][1]]
        else:
            for i in range(1, limit+1):
                primer = sequence[preds[n][0] : preds[n][1]+i*param]
                Tm = primer3.calc_tm(primer)
                if abs(Tm - median) <= Tm_range:
                    break
            Tms[n] = Tm
            preds[n] = [preds[n][0], preds[n][1]+i*param]
    return preds


def assembly_design(sequence, Tm_range, expect_length, expect_overlap, flexibility):
    infor = cutoff(sequence, expect_length-expect_overlap, expect_overlap)
    data = data_processing(infor, flexibility, expect_overlap)
    preds = predict(data)
    try:
        preds = adjust_Tm(sequence, preds, Tm_range)
        return idx_to_oligo(sequence, preds)
    except Exception as ex:
        traceback.print_exc()
        print('Wrong!')
        exit()

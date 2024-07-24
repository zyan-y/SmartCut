

import numpy as np
import tensorflow as tf
from tensorflow.keras.models import load_model

AlphaBet = {'A':0, 'C':1, 'G':2, 'T':3, '#':4}
DefaultLength = 18


def OneHot(seqs):
    data = np.zeros((len(seqs),DefaultLength,len(AlphaBet.keys())))
    for i in range(len(seqs)):
        seq = seqs[i]
        for j in range(len(seq)):
            data[i][j][AlphaBet[seq[j]]] = 1
    return data

def process_data(seq, length=DefaultLength):
    if len(seq) >= length:
        return seq[0:length//2] + seq[(-1)*length//2:]
    else:
        left = (length - len(seq)) // 2
        right = length - len(seq) - left
        return left*'#' + seq + right*'#'

def predict_data(seqs):
    x = OneHot([process_data(i) for i in seqs])
    model = tf.keras.models.load_model(r'./SmartCut/model.h5')
    y_prob = model.predict(x)
    return y_prob

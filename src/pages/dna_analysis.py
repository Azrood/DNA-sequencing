from .dna_func import *
from django.conf import settings
from Bio import SeqIO
from collections import Counter
import plotly.graph_objects as go
import plotly
from Bio.SeqUtils import GC
from io import StringIO
import matplotlib.pyplot as plt
import matplotlib
import numpy as np


def analysis(filename):
     
    with open(filename, 'r') as f:
        dna_file = SeqIO.parse(f,"fasta")
        l = [record for record in dna_file] 
        dna_record = l[0]
        return dna_record


def delta(x,y):
    return 0 if x == y else 1


def M(seq1,seq2,i,j,k):
    return sum(delta(x,y) for x,y in zip(seq1[i:i+k],seq2[j:j+k]))


def makeMatrix(seq1,seq2,k):
    n = len(seq1)
    m = len(seq2)
    return [[M(seq1,seq2,i,j,k) for j in range(m-k+1)] for i in range(n-k+1)]


def plotMatrix(M,t, seq1, seq2, nonblank = chr(0x25A0), blank = ' '):
    print(' |' + seq2)
    print('-'*(2 + len(seq2)))
    for label,row in zip(seq1,M):
        line = ''.join(nonblank if s < t else blank for s in row)
        print(label + '|' + line)




def analysis_plot(seq):
    dna_freq = Counter(seq)
    fig = go.Figure(
    data=[go.Bar(x=list(dna_freq.keys()),y=list(dna_freq.values()))],
    layout_title_text="Graphe de fréquence"
    )
    return fig

def dotplotx(seq1,seq2):
    fig =plt.figure()
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    x=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # # on y-axis list all sequences of seq 1
    y=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    imgdata = StringIO()
    fig.savefig(imgdata, format='svg')
    imgdata.seek(0)
    data = imgdata.getvalue()
    return data
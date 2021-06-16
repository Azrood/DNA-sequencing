from .dna_func import *
from django.conf import settings
from Bio import SeqIO
from collections import Counter
import plotly.graph_objects as go
import plotly
from Bio.SeqUtils import GC


def analysis(filename):
     
    file = open(filename, 'r')
    dna_file = SeqIO.parse(file,"fasta")
    l = [record for record in dna_file] 
    dna_record = l[0]
    return dna_record

def analysis_plot(seq):
    dna_freq = Counter(seq)
    fig = go.Figure(
    data=[go.Bar(x=list(dna_freq.keys()),y=list(dna_freq.values()))],
    layout_title_text="Graphe de fréquence"
    )
    return fig

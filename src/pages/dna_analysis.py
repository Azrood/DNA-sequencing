from .dna_func import *
from django.conf import settings
from Bio import SeqIO
from collections import Counter
import plotly.graph_objects as go
import plotly

file = open(f'{settings.MEDIA_ROOT}/covid19.fasta','r')
dna_file = SeqIO.parse(file,"fasta")
l = [record for record in dna_file] 

def analysis():
     
    dna_record = l[0]
    return dna_record

def analysis_plot():
    dna_freq = Counter(analysis().seq)
    fig = go.Figure(
    data=[go.Bar(x=list(dna_freq.keys()),y=list(dna_freq.values()))],
    layout_title_text="test"
    )
    return fig

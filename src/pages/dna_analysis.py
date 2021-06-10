from .dna_func import *
from django.conf import settings
from Bio import SeqIO

def analysis():
    file = open(f'{settings.MEDIA_ROOT}/covid19.fasta','r')
    dna_file = SeqIO.parse(file,"fasta")
    #st.write(dna_record)
    l = [record for record in dna_file]  
    dna_record = l[0]
    dna_seq= dna_record.seq


    return dna_record.description
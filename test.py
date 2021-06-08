import streamlit as st
from Bio.Seq import Seq
from Bio import SeqIO
import neatbio.sequtils as utils
from collections import Counter
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import numpy as np
import streamlit.components.v1 as components
from stmol import component_3dmol
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import pandas as pd

LOGGED = False

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


def dotplot(seq1,seq2,k = 1,t = 1):
    M = makeMatrix(seq1,seq2,k)
    plotMatrix(M, t, seq1,seq2) #experiment with character choice


def dotplotx(seq1,seq2):
    fig,ax = plt.subplots()
    
    ax.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.rcParams.update({'font.size': 5})
    plt.show()
    return fig



def main():
    """Bioinformatics Genome analysis web app"""
    st.title("DNA Genome analysis and Cosine Similarity Analysis web application")
    menu= ["Introduction","DNA sequence Analysis","Dotplot Analysis","3D Visualization","About us"]

    choice= st.sidebar.selectbox("Select Option",menu)

    if choice=="Introduction":
        st.subheader("Welcome to our Sequence Analysis Application :)")

    elif choice=="DNA sequence Analysis":
        
        st.subheader("DNA sequence Analysis will be done here.")
        seq_file=st.file_uploader("Upload the .FASTA file for any DNA analysis of the considered Genome.", type=["fasta","fa"])
                
        if seq_file is not None:
            with open('test.fasta',"wb") as f: 
                for line in seq_file.readlines():
                    f.write(line)
            file = open('test.fasta','r')
            dna_file = SeqIO.parse(file,"fasta")
            #st.write(dna_record)
            l = [record for record in dna_file]  
            dna_record = l[0]
            dna_seq= dna_record.seq
            cols=st.beta_columns(2)
            details= cols[0].radio("Details of the DNA as provided by NCBI database:",("DNA Record description", "Sequence"))
            if details=="DNA Record description":
                cols[0].write(dna_record.description)
            elif details=="Sequence":
                cols[0].code(dna_record.seq)


            #Nucleotide
            cols[1].subheader("Nucleotide Frequency :")
            dna_freq=Counter(dna_seq)
            cols[1].write(dna_freq)
            tpCols=st.beta_columns(4)
            adenine_color=tpCols[0].color_picker("Toggle the Adenine Colour ")
            guanine_color=tpCols[1].color_picker("Toggle the Guanine Colour ")
            thymine_color=tpCols[2].color_picker("Toggle the Thymine Colour ")
            cytosine_color=tpCols[3].color_picker("Toggle the Cytosine Colour ")


            if st.button("Plot frequency"):
                fig,ax = plt.subplots()
                barlist=ax.bar(dna_freq.keys(),dna_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(guanine_color)
                barlist[2].set_color(thymine_color)
                barlist[3].set_color(cytosine_color)
                st.pyplot(fig)


            st.subheader("DNA complete Composition")

            gc_score= utils.gc_content(str(dna_seq))
            at_score=utils.at_content(str(dna_seq))
            st.json({"GC Content(for heat stability)": gc_score,"AT Content":at_score })

            #protein synthesis
            st.subheader("Protein Synthesis operations on the DNA :")
            p1=dna_seq.translate()
            aa_freq= Counter(str(p1))
            if st.checkbox("Transcription :"):
                st.write(dna_seq.transcribe())
            elif st.checkbox("Translation :"):
                st.write(dna_seq.translate())
            elif st.checkbox("Complement :"):
                st.write(dna_seq.complement())
            elif st.checkbox("Amino Acid frequency :"):
                st.write(aa_freq)

            elif st.checkbox("Plot the Amino Acid frequency :"):
                aa_color=st.color_picker("Pick the Amino acid color:")
                #barlist= plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                #barlist[2].set_color(aa_color)
                f,a =plt.subplots()
                a.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot(f)

            elif st.checkbox("The complete Amino acid name is given as"):
                amino_acids = [a for a in p1.split('*')]
                amino = [utils.get_acid_name(utils.convert_1to3(a)) for a in amino_acids ]
                df = pd.DataFrame({'amino_acids':amino_acids,'full':amino})
                df['count'] = df['amino_acids'].apply(len)
                st.write(df.nlargest(20,'count'))
                # aa_name= str(p1).replace("*","")
                # aa3= utils.convert_1to3(aa_name)
                # st.write(aa_name)
                # st.write("========================")
                # st.write(aa3)
                # st.write("========================")
                # st.write(utils.get_acid_name(aa3))



    elif choice=="Dotplot Analysis":
        st.subheader("Generate Dotplot for the comparision between two DNA sequences here.")
        try:
            seq_file1=st.file_uploader("Upload the first .FASTA file for any DNA analysis of the considered Genome.", type=["fasta","fa"])
            with open('test1.fasta',"wb") as f: 
                for line in seq_file1.readlines():
                    f.write(line)
            seq_file2=st.file_uploader("Upload the second .FASTA file for any DNA analysis of the considered Genome.", type=["fasta","fa"])

            with open('test2.fasta',"wb") as f: 
                for line in seq_file2.readlines():
                    f.write(line) 
        except:
            pass


        if seq_file1 and seq_file2 is not None:

            file1 = open('test1.fasta','r')
            file2 = open('test2.fasta','r')
            dna_record1= list(SeqIO.parse(file1,"fasta"))
            dna_record2= list(SeqIO.parse(file2,"fasta"))
            # seq1 = Seq('ATGTCGATGACGCT')
            # seq2 = Seq('GTACTGCAGTC')
            #st.write(dna_record)
            choice1 = st.selectbox("Selectionnez votre sequence",[rec for rec in dna_record1],format_func=lambda x: str(x.description))
            choice2 = st.selectbox("Selectionnez votre sequence",[rec for rec in dna_record2],format_func=lambda x: str(x.description))
            dna_seq1= choice1.seq
            dna_seq2= choice2.seq
            # print(dna_seq1)
            # print(dna_seq2)
            


            details= st.radio("Details of the DNA as provided by NCBI database:",("Record details from the NCBI database", "Gene Sequence"))
            if details=="Record details from the NCBI database":
                st.code(choice1.description)
                st.write("===And the other Record is decribed as :===")
                st.code(choice1.description)

            elif details=="Gene Sequence":
                st.write(choice1.seq)
                st.write("===And the other sequence can be given as: ===")
                st.write(choice1.seq)
            st.subheader('Alignement des sequences')
            alignments = pairwise2.align.globalxx(dna_seq1[0:100],dna_seq2[0:100])
            st.code (format_alignment(*alignments[0]))

            st.subheader('Dotplot')

            display_limit=st.number_input("Select maximum number of Nucleotides",10,1000,50)
            if st.button("Push here for Dotplot :)"):
                st.write("Comparing the first {} nucleotide of the two sequences".format(display_limit))
                t=dotplotx(dna_seq1[0:display_limit],dna_seq2[0:display_limit])
                st.pyplot(t)
                # try:
                #     st.pyplot()
                # except:
                #     pass

    # if LOGGED :
    #     block1=st.empty()
    #     block2=st.empty()
    #     block3=st.empty()
        
    #     username = block1.text_input('Username')
    #     password = block2.text_input('Password')
    #     if block3.button('Login'):
    #         if True:
    #             block1.empty()
    #             block2.empty()
    #             block3.empty()

    #             deploy()
    # else:
    #     deploy()

if __name__=='__main__':
    main()
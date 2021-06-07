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
    plt.imshow(np.array(makeMatrix(seq1,seq2,1)))
    # on x-axis list all sequences of seq 2
    xt=plt.xticks(np.arange(len(list(seq2))),list(seq2))
    # on y-axis list all sequences of seq 1
    yt=plt.yticks(np.arange(len(list(seq1))),list(seq1))
    plt.show()



def main():
    """Bioinformatics Genome analysis web app"""


    st.title("DNA Genome analysis and Cosine Similarity Analysis web application")
    menu= ["Introduction","DNA sequence Analysis","Dotplot Analysis","About us"]
    choice= st.sidebar.selectbox("Select Option",menu)

    if choice=="Introduction":
        st.subheader("Welcome to our Sequence Analysis Application :)")
        components.html('''TEST''')
        component_3dmol()
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

            details= st.radio("Details of the DNA as provided by NCBI database:",("DNA Record description", "Sequence"))
            if details=="DNA Record description":
                st.write(dna_record.description)
            elif details=="Sequence":
                st.write(dna_record.seq)

            #Nucleotide
            st.subheader("Nucleotide Frequency :")
            dna_freq=Counter(dna_seq)
            st.write(dna_freq)
            adenine_color=st.color_picker("Toggle the Adenine Colour ")
            guanine_color=st.color_picker("Toggle the Guanine Colour ")
            thymine_color=st.color_picker("Toggle the Thymine Colour ")
            cytosine_color=st.color_picker("Toggle the Cytosine Colour ")


            if st.button("Plot frequency"):
                barlist=plt.bar(dna_freq.keys(),dna_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(guanine_color)
                barlist[2].set_color(thymine_color)
                barlist[3].set_color(cytosine_color)
                st.pyplot()


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
                plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot()

            elif st.checkbox("The complete Amino acid name is given as"):
                aa_name= str(p1).replace("*","")
                aa3= utils.convert_1to3(aa_name)
                st.write(aa_name)
                st.write("========================")
                st.write(aa3)



                st.write("========================")
                st.write(utils.get_acid_name(aa3))



    elif choice=="Dotplot Analysis":
        st.subheader("Generate Dotplot for the comparision between two DNA sequences here.")
        seq_file1=st.file_uploader("Upload the first .FASTA file for any DNA analysis of the considered Genome.", type=["fasta","fa"])
        seq_file2=st.file_uploader("Upload the second .FASTA file for any DNA analysis of the considered Genome.", type=["fasta","fa"])


        if seq_file1 and seq_file2 is not None:
            with open('test1.fasta',"wb") as f: 
                for line in seq_file1.readlines():
                    f.write(line)
            with open('test2.fasta',"wb") as f: 
                for line in seq_file2.readlines():
                    f.write(line)    

            file1 = open('test1.fasta','r')
            file2 = open('test2.fasta','r')
            dna_record1= SeqIO.read(file1,"fasta")
            dna_record2= SeqIO.read(file2,"fasta")
            #st.write(dna_record)
            dna_seq1= dna_record1.seq
            dna_seq2= dna_record2.seq

            details= st.radio("Details of the DNA as provided by NCBI database:",("Record details from the NCBI database", "Gene Sequence"))
            if details=="Record details from the NCBI database":
                st.write(dna_record1.description)
                st.write("===And the other Record is decribed as :===")
                st.write(dna_record2.description)

            elif details=="Gene Sequence":
                st.write(dna_record1.seq)
                st.write("===And the other sequence can be given as: ===")
                st.write(dna_record2.seq)


            display_limit=st.number_input("Select maximum number of Nucleotides",10,200,50)
            if st.button("Push here for Dotplot :)"):
                st.write("Comparing the first {} nucleotide of the two sequences".format(display_limit))
                dotplotx(dna_seq1[0:display_limit],dna_seq2[0:display_limit])
                st.pyplot()



    elif choice=="About us":
        st.subheader("About the application and about us :)")

if __name__=='__main__':
    main()
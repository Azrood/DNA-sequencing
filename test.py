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


    st.title("Application web sécurisée d'analyse de séquence ADN")
    menu= ["Introduction","Analyse de séquence ADN","Comparaison de séquences ADN","Visualisation 3D","À propos"]
    choice= st.sidebar.selectbox("Sélectionner option",menu)

    if choice=="Introduction":
        st.subheader("Bienvenue chez notre application d'analyse de séquence ADN :)")

    elif choice=="Analyse de séquence ADN":
        st.subheader("Analyse de séquence ADN sera faite ici :")
        seq_file=st.file_uploader("Uploadez le fichier .FASTA pour toute analyse ADN du génome considéré.", type=["fasta","fa"])
                 
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
            details= cols[0].radio("Détails de l'ADN fournit par la base de donnée NCBI:",("Description des enregistrements ADN", "Séquence"))
            if details=="Description des enregistrements ADN":
                cols[0].write(dna_record.description)
            elif details=="Séquence":
                cols[0].code(dna_record.seq)


            #Nucleotide
            cols[1].subheader("Fréquence de nucléotide :")
            dna_freq=Counter(dna_seq)
            cols[1].write(dna_freq)
            tpCols=st.beta_columns(4)
            adenine_color=tpCols[0].color_picker("Activer la couleur Adenine")
            guanine_color=tpCols[1].color_picker("Activer la couleur Guanine")
            thymine_color=tpCols[2].color_picker("Activer la couleur Thymine")
            cytosine_color=tpCols[3].color_picker("Activer la couleur Cytosine")


            if st.button("Graphe de fréquence"):
                fig,ax = plt.subplots()
                barlist=ax.bar(dna_freq.keys(),dna_freq.values())
                barlist[0].set_color(adenine_color)
                barlist[1].set_color(guanine_color)
                barlist[2].set_color(thymine_color)
                barlist[3].set_color(cytosine_color)
                st.pyplot(fig)


            st.subheader("Composition complète de l'ADN")

            gc_score= utils.gc_content(str(dna_seq))
            at_score=utils.at_content(str(dna_seq))
            st.json({"Taux de GC(Stabilité thermique)": gc_score,"Taux de AT":at_score })

            #protein synthesis
            st.subheader("Opérations de synthèses des protéines sur l'ADN:")
            p1=dna_seq.translate()
            aa_freq= Counter(str(p1))
            if st.checkbox("Transcription :"):
                st.write(dna_seq.transcribe())
            elif st.checkbox("Traduction :"):
                st.write(dna_seq.translate())
            elif st.checkbox("Complément :"):
                st.write(dna_seq.complement())
            elif st.checkbox("Fréquence des acides aminées :"):
                st.write(aa_freq)

            elif st.checkbox("Graphe de fréquence des acides aminées :"):
                aa_color=st.color_picker("Choisir la couleur des acides aminées:")
                #barlist= plt.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                #barlist[2].set_color(aa_color)
                f,a =plt.subplots()
                a.bar(aa_freq.keys(),aa_freq.values(),color=aa_color)
                st.pyplot(f)

            elif st.checkbox("Le nom complet de l'acide aminée est :"):
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



    elif choice=="Comparaison de séquences ADN":
        st.subheader("Génère une comparaison entre 2 séquences ADN..")
        try:
            seq_file1=st.file_uploader("Uploader le premier fichier .FASTA pour toute analyse ADN du génome considéré.", type=["fasta","fa"])

            with open('test1.fasta',"wb") as f: 
                for line in seq_file1.readlines():
                    f.write(line)
            seq_file2=st.file_uploader("Uploader le second fichier .FASTA pour toute analyse ADN du génome considéré.", type=["fasta","fa"])

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
            


            details= st.radio("Détails de l'ADN fournis par la base de donnée NCBI :",("Détails des enregistrement de la base de donnée NCBI", "Séquence de gènes"))
            if details=="Détails des enregistrement de la base de donnée NCBI":
                st.code(choice1.description)
                st.write("===L'autre enregistrement peut être décrit comme :===")
                st.code(choice2.description)

            elif details=="Séquence de gènes":
                st.write(choice1.seq)
                st.write("===L'autre séquence peut être décrite comme: ===")
                st.write(choice2.seq)
            st.subheader('Alignement des séquences')
            alignments = pairwise2.align.globalxx(dna_seq1[0:100],dna_seq2[0:100])
            st.code (format_alignment(*alignments[0]))

            st.subheader('Dotplot')

            display_limit=st.number_input("Selectionnez le nombre maximum de nucléotides",10,200,50)
            if st.button("Cliquez ici pour le graphe de comapraison :)"):
                st.write("Comparaison des premières {} nucléotide des 2 séquences".format(display_limit))
                t=dotplotx(dna_seq1[0:display_limit],dna_seq2[0:display_limit])
                st.pyplot(t)


    elif choice=="Visualisation 3D":
        component_3dmol()
    elif choice=="À propos":
        st.subheader("About the application and about us :)")
        st.write("Application réalisée par : CHEIKH Mohammed Nabil et L'HICHOU Anas. Sous la supervision de Pr. BERQIA")

if __name__=='__main__':
    main()
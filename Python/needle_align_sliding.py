"""needle_align_sliding.py: Script to pairwise align sequence read and generate per-nucleotide scores using a sliding window."""
__author__ = "Devang Mehta"
__copyright__ = "Copyright 2017, ETH Zurich"
__version__ = "1.0.0"
__maintainer__ = "Devang Mehta"
__email__ = "devang@ethz.ch"

import sys
import os
import glob
import click
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import align, format_alignment
from Bio.SubsMat import MatrixInfo as matlist
from Bio.Alphabet import IUPAC

@click.command()
@click.argument('inputdir', type=click.Path(exists=True, readable=True))
@click.argument('refseq', type=click.Path(exists=True, readable=True))

def main(inputdir,refseq,**kwargs):
    target=inputdir+"/result/"
    os.makedirs(target, exist_ok=True)
    ref=SeqIO.read(refseq,"fasta")
    refs=ref.seq
    print("-------------------------------------------------Reading reference \n")
    files=glob.glob(inputdir+"/*.fasta")
    for ffile in files:
        print(os.path.splitext(os.path.basename(ffile))[0],"\n")
        fname=os.path.abspath(inputdir)+"/result/"+os.path.splitext(os.path.basename(ffile))[0]+".csv"
        fout=open(fname,"w")
        num=0
        fout.write("height"+","+"sequence")
        nuc=list(range(1,2692))
        mylist=','.join(map(str, nuc))
        fout.write(","+mylist)
        for record in SeqIO.parse(ffile, 'fasta'):
            num+=1
            print("--------Reading record \n")
            records=record.seq
            gap_open = -10
            gap_extend = -0.5
            print("--------Aligning \n")
            alns = pairwise2.align.globalxs(records,refs,gap_open,gap_extend)
            best_aln = alns[0]
            aligned_A, aligned_B, score, begin, end = best_aln
            sa, sb, sl = aligned_A, aligned_B, len(aligned_A)
            fout.write('\n'+str(num)+","+record.id)
            for i in range (5,sl-6):
                sa1=sa[i-5:i+6]
                print(sa1)
                sb1=sb[i-5:i+6]
                print(sb1+'\n')
                matches=float(0)
                for ii in range(11):
                    if sa1[ii] == sb1[ii]:
                        matches+=1
                seq_id=100-((100*matches)/ 11)
                print(matches)
                print(seq_id)
                fout.write(','+str(seq_id))

            print("-----------------------------------------------------------------------------------------------\n\n\n")
            fout.write("\n")


        fout.close()

if __name__ == '__main__':

	main()

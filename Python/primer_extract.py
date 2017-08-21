"""primer_extract.py: Script to extract a sequence lying between user-defined primer sequences."""
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

@click.command()
@click.argument('inputdir', type=click.Path(exists=True, readable=True))
@click.argument('primer_f',  type=click.STRING)
@click.argument('primer_r',  type=click.STRING)
#100for: GTTTATTTC
#100rev: ACAGACAAG

#wholefor: CCTTTAATTTGAAC
#wholerev: TTGCTTTTCC


def main(inputdir,primer_f, primer_r):

    target=inputdir+"/output/"
    os.makedirs(target, exist_ok=True)

    files=glob.glob(inputdir+"/*.fasta")
    for ffile in files:
        fname=os.path.abspath(inputdir)+"/output/"+os.path.splitext(os.path.basename(ffile))[0]+".fasta"
        fout=open(fname,"w")
        for record in SeqIO.parse(ffile, 'fasta'):
            finder(fout,record.id,record.seq,primer_f,primer_r)
            finder(fout,record.id+"_rev",record.seq.reverse_complement(),primer_f,primer_r)

        fout.close()


def finder(fout,id,seq,primer_f,primer_r):
    spos=0
    while spos >=0:
        spos=seq.find(primer_f,spos)
        if spos >=0:
            tpos=seq.find(primer_r,spos)
            if tpos >=0:
                print('>'+id+" | "+str(spos+len(primer_f))+":"+str(tpos))
                fout.write('>'+id+" | "+str(spos+len(primer_f))+":"+str(tpos)+"\n")
                fout.write(str(seq[(spos+len(primer_f)):tpos])+"\n")
                spos+=1
            else:
                spos=-1



if __name__ == '__main__':

	main()

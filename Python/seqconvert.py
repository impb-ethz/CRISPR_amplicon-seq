"""sequconvert.py: Script to convert sequences from format to format"""
__author__ = "Devang Mehta"
__copyright__ = "Copyright 2017, ETH Zurich"
__version__ = "1.0.0"
__maintainer__ = "Devang Mehta"
__email__ = "devang@ethz.ch"
import sys
import glob
import os
import click
from Bio import SeqIO

@click.command()
@click.argument('inputdir', type=click.Path(exists=True, readable=True))
@click.option('--informat', default='fastq',type=click.Choice(['fasta','fastq','tab','gb']), help='Input-file format (default=fastq).')
@click.option('--outformat',default='fasta',type=click.Choice(['fasta','tab','gb', 'fastq','txt']), help='Input-file format (default=fasta).')

def main(inputdir,informat,outformat):

    target=inputdir+"/output/"+outformat
    os.makedirs(target, exist_ok=True)

    files=glob.glob(inputdir+"/*."+informat)
    for ffile in files:
        fname=os.path.abspath(inputdir)+"/output/"+outformat+"/"+os.path.splitext(os.path.basename(ffile))[0]+"."+outformat
        fout=open(fname,"w") #overwrite
        print(inputdir)
        print(os.path.splitext(os.path.basename(ffile))[0])
        print(fname)
        #read/write fasta file
        for record in SeqIO.parse(ffile,informat):
            SeqIO.write(record,fout,outformat)
        fout.close()

if __name__ == '__main__':
    main()

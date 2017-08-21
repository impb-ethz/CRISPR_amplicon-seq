"""ext_replace.py: Script to change file extensions."""
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

def main(inputdir,**kwargs):

    files=glob.glob(inputdir+"/*.fasta")
    for ffile in files:
        base=os.path.splitext(ffile)[0]
        print(base)
        os.rename(ffile, base+".csv")

if __name__ == '__main__':
    main()

#!/usr/bin/env python
"""annotate.py: Script to 

1. blast dna segments against a reference using tblastn
2. translate hsps into protein
3. detect MATCH, SNP and INDEL or resulting protein using Muscle Alignment

"""
__author__ = "Matthias Hirsch-Hoffmann, Devang Mehta"
__copyright__ = "Copyright 2017, ETH Zuerich"
__version__ = "1.0.0"
__maintainer__ = "Matthias Hirsch-Hoffmann"
__email__ = "hirschhm@ethz.ch"
__status__ = "Production"

import os,sys
import uuid
import tempfile
import json
import click
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML
from Bio.Align.Applications import MuscleCommandline
#define excutables for blast and muscle
blast_exe = 'tblastn'
muscle_exe = 'muscle3.8.31_i86linux64'

@click.command()
#arguments=mandatory, option=possible
#*************** ARGUMENETS *****************************
#input-file in fasta format
@click.argument('inputfile', type=click.Path(exists=True,readable=True))
#blastdatabase
# must exists in fasta format and as formated database
@click.argument('blastdatabase', type=click.Path(exists=True,readable=True))


def annotate(inputfile,blastdatabase):
	#stores information from the blastdb
	protein={}
	#read blastdb fasta file
	for record in SeqIO.parse(blastdatabase, 'fasta'):
		#create new result entry
		protein[record.id]=str(record.seq).replace('*','')
		#cleanup result directroy/file
		if os.path.isfile(inputfile+"."+record.id+".txt"):
			os.remove(inputfile+"."+record.id+".txt")	
	
	#create temporary input/output directory
	tmpdir=tempfile.mkdtemp()
	#this will hold the results	
	rdata={}
	#read file and annotate
	for record in SeqIO.parse(inputfile, 'fasta'):
		#create new result entry
		rdata[record.id]={}
		rdata[record.id]['seq']=str(record.seq)
		rdata[record.id]['length']=len(record.seq)
		rdata[record.id]['proteins']={}
		#print(record.id,'-',len(record.seq))
	
		#define blast base file
		blastfile=str(tmpdir)+"/"+str(uuid.uuid4())
	
		#write blast input file
		with open(blastfile+'.in', "wt") as output_handle:
			SeqIO.write(record,output_handle,'fasta')
		output_handle.close()
		#build tblastn command
		cline=NcbitblastnCommandline(
				    query = blastdatabase
				, subject = blastfile+".in" 
				,  evalue = 0.01
				,     out = blastfile+".out"
				,     cmd = blast_exe
				,    task = 'tblastn'
				,  outfmt = 5) #xml
		#execute command
		stout, stderr = cline()
		#dictionary for results
		aligns={}
		#read result
		for b_record in NCBIXML.parse(open(blastfile+'.out', 'r')):
			#read all alignments
			for alignment in b_record.alignments:
				#read all hsp 
				for hsp in alignment.hsps:
					#query-frame and subject-frame
					qframe,sframe = hsp.frame
					#does align key exists in aligns?
					if b_record.query in aligns.keys():
						#append alignment
						aligns[b_record.query].append([sframe,hsp.sbjct_start,hsp.sbjct_end,hsp.query_start,hsp.query_end,hsp.expect])
					else:
						#create key
						aligns[b_record.query]={}
						#all alignment
						aligns[b_record.query]=[[sframe,hsp.sbjct_start,hsp.sbjct_end,hsp.query_start,hsp.query_end,hsp.expect]]
	
		#remove input/output files
		os.remove(blastfile+'.in')
		os.remove(blastfile+'.out')
					
		#evaluate results
		#loop through all keys
		for key in sorted(aligns):
			rdata[record.id]['proteins'][key]={}
			#print(key,len(aligns[key]),aligns[key])
			rdata[record.id]['proteins'][key]['hsps']=aligns[key]
			#def and set result variables
			minprotpos=len(record.seq)+99
			maxprotpos=1
			minnucpos=-1
			maxnucpos=-1
			strand=0
			#loop through aligns
			for i in range(0,len(aligns[key])):
				#extract elements
				f2_frame,f2_start,f2_end,q2_start,q2_end,evalue=aligns[key][i]
				
				#check same frame direction
				if f2_frame < 0 and strand==0:
					#reverse strand and strand not yet set
					strand=-1
				elif f2_frame > 0 and strand==0:
					#forward strand and strand not yet set
					strand=1
				elif f2_frame > 0 and strand>0:
					#ok - forward strand
					nop=0
				elif f2_frame < 0 and strand<0:
					#ok - reverse strand
					nop=0
				else:
					#strand contradicting, break here, the sequence will not be further processed
					#print("strand error")
					strand=-99
					break
				#reverse, we have to switch start and end 
				if f2_frame < 0:
					#switch nuc start-end
					f2_start,f2_end=f2_end,f2_start
				
				#def smallest protein position
				if q2_start < minprotpos:
					minprotpos=q2_start
					#with that, set smalles nucleotide position
					minnucpos=f2_start
				#def largest protein position
				if q2_end > maxprotpos:
					maxprotpos=q2_end
					#with that, set largest nucleotide position
					maxnucpos=f2_end
	
			#strand and protein start-end position are ready
			rdata[record.id]['proteins'][key]['strand']=strand
			rdata[record.id]['proteins'][key]['minprotpos']=minprotpos
			rdata[record.id]['proteins'][key]['maxprotpos']=maxprotpos
	
			#evaluate 
			if strand < 0:
				#reverse minnucposition in result is current maxnucpos
				rdata[record.id]['proteins'][key]['minnucpos']=maxnucpos
				#check if we run over breakpoint
				if minnucpos > maxnucpos:
					#normal piece, just turn minnucpos and maxnucpos
					#print("ok reverse")
					#print(strand,maxnucpos,minnucpos,minprotpos,maxprotpos)		
					rdata[record.id]['proteins'][key]['case']='ok reverse'
					rdata[record.id]['proteins'][key]['maxnucpos']=minnucpos
				else:
					#over breakpoint
					print('case not handled')
					sys.exit(1)
					#print("over breakpoint reverse")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					#print(strand,maxnucpos,(len(record.seq)+minnucpos),minprotpos,maxprotpos)
					rdata[record.id]['proteins'][key]['case']='over breakpoint reverse'
					rdata[record.id]['proteins'][key]['maxnucpos']=(len(record.seq)+minnucpos)
			else:
				rdata[record.id]['proteins'][key]['minnucpos']=minnucpos
				#check if we run over breakpoint
				if minnucpos > maxnucpos:
					#over breakpoint
					print('case not handled')
					sys.exit(1)
					#print("over breakpoint forward")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					#print(strand,minnucpos,(len(record.seq)+maxnucpos),minprotpos,maxprotpos)
					rdata[record.id]['proteins'][key]['case']='over breakpoint forward'
					rdata[record.id]['proteins'][key]['maxnucpos']=(len(record.seq)+maxnucpos)
				else:
					#normal piece
					#print("ok forward")
					#print(strand,minnucpos,maxnucpos,minprotpos,maxprotpos)		
					rdata[record.id]['proteins'][key]['case']='ok forward'
					rdata[record.id]['proteins'][key]['maxnucpos']=maxnucpos

			sequence={}

			for (fr,sp,ep,ps,pe,sc) in rdata[record.id]['proteins'][key]['hsps']:
				partseq=rdata[record.id]['seq'][sp:(ep)]
				#print(partseq)
				if strand < 0:
					tmp=Seq(partseq)
					partseq=str(tmp.reverse_complement())	
				#print(partseq)
				partseq=seqN(partseq)
				#print(partseq)
				#partseq must me a multiple of 3
				sequence[ps]=partseq
			
			transseq=''
			for i in sorted(sequence):
				transseq+=sequence[i]
				

			#tmp=Seq(rdata[record.id]['proteins'][key]['cds'])
			tmp=Seq(transseq)
			rdata[record.id]['proteins'][key]['pep']=str(tmp.translate()).replace('X','')
#analysis of results
	
	for id in rdata:
		for p in rdata[id]['proteins']:
#"strand": -1, "minprotpos": 1, "maxprotpos": 136, "minnucpos": 2254, "case": "ok reverse", "maxnucpos": 2647}
			outputfile=inputfile+"."+p+".txt"
			fout= open(outputfile,'at')
			
			alignres=seqMuscle(tmpdir,rdata[id]['proteins'][p]['pep'],protein[p])
			#print('read:'+alignres['read'])
			#print('ref :'+alignres['ref'])
			
			if alignres['read']==alignres['ref']:
				#print(id+"\tMATCH")
				fout.write(id+"\tMATCH")
			elif len(alignres['read'])==len(alignres['ref']) and alignres['read'].find('-')==alignres['ref'].find('-'):
				#print(id+"\tSNP")
				fout.write(id+"\tSNP")
				snp=defSNP(alignres['read'],alignres['ref'])
				#print(snp)
				#print(len(snp))
				fout.write("\t\t"+str(len(snp)))
				for s in snp:
					#print(s)
					fout.write("\t"+s)
			else:
				#print(id+"\tINDEL")
				fout.write(id+"\tINDEL")
				gappos=alignres['read'].find('-')
#				print(gappos)
				refgappos=alignres['ref'].find('-')
#				print(refgappos)
				if gappos > -1:
					if refgappos > -1 and refgappos < gappos:
						gappos=refgappos
				elif refgappos > -1:
					gappos = refgappos
				
				#print(gappos)			
				fout.write("\t"+str(gappos))
				if gappos > 0:
					#print(alignres['read'][0:gappos])
					#print(alignres['ref'][0:gappos])
					#print(gappos)
					snp=defSNP(alignres['read'][0:gappos],alignres['ref'][0:gappos])
					#print(snp)
					#print(len(snp))
					fout.write("\t"+str(len(snp)))
					for s in snp:
						#print(s)
						fout.write("\t"+s)
			
			fout.write('\n')
			fout.close()
			#print(alignres)
			
	#write result
#	outputfile="./"+os.path.splitext(os.path.basename(inputfile))[0]
#	#print(json.dumps(rdata))
#	with open(outputfile+".json",'w') as fout:
#		json.dump(rdata,fout)
#	fout.close()

	os.rmdir(tmpdir)
	#return genome
	return outputfile

def defSNP(readseq,refseq):
	snp=[]
	for i in range(0,len(readseq)):
		#print(i)
		read=readseq[i:i+1]
		ref=refseq[i:i+1]
		#print(read+str(i)+ref)
		if ref!=read:
			snp.append(ref+str(i)+read)	
	return snp	

def seqMuscle(tmpdir,readseq,refseq):
	musclefile = tmpdir+"/"+str(uuid.uuid4())
	#define in 
	in_file = musclefile+".in"
	# define out file
	out_file = musclefile+".out"

	# write in_file fasta file
	fout=open(in_file,'wt')
	fout.write('>read\n') #submitted record id
	fout.write(readseq+'\n') #submitted record sequence
	fout.write('>ref\n') #target genome id
	fout.write(refseq+'\n') #target genome sequence * 3
	fout.close()
	#execute muscle
	muscle_cline = MuscleCommandline(muscle_exe, input=in_file, out=out_file)
	stdout, stderr = muscle_cline()	
	#read result
	resseq={} #fill result into dict to have direct access to seqrecord
	#read out_file file
	for resrecord in SeqIO.parse(out_file, "fasta"):
		#print(resrecord.seq)
		resseq[str(resrecord.id)]=str(resrecord.seq)
	#delete temporary input/output files
	os.remove(in_file)
	os.remove(out_file)
	#return alignment
	return resseq

	
def seqN(sequence):
	#fills sequence with N to be a multiple of 3
	mod = len(sequence)%3
	if mod==1:
		sequence+='NN'
	elif mod==2:
		sequence+='N'
	return sequence	
				
if __name__ == '__main__':
	annotate()

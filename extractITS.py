#!/usr/bin/env python

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import argparse
import subprocess
import sys
import pandas as pd
import os

def get_params(argv):
	parser = argparse.ArgumentParser(description='Extract ITS1 from fungal genome rapidly')
	parser.add_argument('-i', '--i', help="input genome file", required=True)
	parser.add_argument('-o', '--o', help="output directory", required=True)
	parser.add_argument('-which', '--which', help="Which ITS sequence to extract (ITS1|ITS2) default=ITS1", default='ITS1')
	parser.add_argument('-cpu', '--cpu', help="number of threads/cores to use", required=False, default='48')
	parser.add_argument('-name', '--name', help="name", required=False, default='genome')
	a = parser.parse_args()
	return a

def terminal(cmd):
	p = subprocess.Popen(cmd, shell=True, stdout = subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
	p.wait()



if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	if a.o[-1]!='/':
		a.o=a.o+'/'

	if not os.path.exists(a.o):
		os.makdirs(a.o)

	if a.name=='genome':
		name=a.i.split('/')[-1].split('.')[0]
	else:
		name=a.name

	if a.which not in ['ITS1','ITS2']:
		print("ERROR, argument --which not recognized. Please choose between 'ITS1' and 'ITS2'")

	## Running barrnap
	terminal('barrnap --kingdom euk --threads '+a.cpu+' '+a.i+' > '+a.o+'barrnap.tmp')

	## Parsing and getting rDNA gene cluster coordinates
	gff=pd.read_csv(a.o+'barrnap.tmp', header=None,comment='#',sep='\t')
	gff=gff[gff[8].str.contains('18S') | gff[8].str.contains('28S')].sort_values(by=[0,6,3,4]).reset_index(drop=True)
	regions=[]
	print(gff)
	for i in gff.index:
		try:
			if 'partial' in gff.loc[i,8]:
				pass
			elif gff.loc[i,6]=='+':
				if '18S_rRNA' in gff.loc[i,8] and '28S_rRNA' in gff.loc[i+1,8] and gff.loc[i,0]==gff.loc[i+1,0] and gff.loc[i,6]==gff.loc[i+1,6]:
					regions.append((gff.loc[i,0],gff.loc[i,6],gff.loc[i,3]-1,gff.loc[i+1,4]-1))
			elif gff.loc[i,6]=='-':
				if '28S_rRNA' in gff.loc[i,8] and '18S_rRNA' in gff.loc[i+1,8] and gff.loc[i,0]==gff.loc[i+1,0] and gff.loc[i,6]==gff.loc[i+1,6]:
					regions.append((gff.loc[i,0],gff.loc[i,6],gff.loc[i,3]-1,gff.loc[i+1,4]-1))
			else:
				print('Error, impossible to read '+gff.loc[i,6])
		except:
			pass

	#terminal('rm '+a.o+'barrnap.tmp')
	print(regions)

	## Opening genome fasta, extracting regions
	with open(a.i,'r') as handle:
		fasta=SeqIO.to_dict(SeqIO.parse(handle,'fasta'))

	extractedRegions=[]
	for r in regions:
		print(r)
		SEQ=fasta[r[0]][r[2]:r[3]]
		SEQ.id=r[0]+'_'+r[1]+'_'+str(r[2])+'_'+str(r[3])
		SEQ.name=SEQ.id
		SEQ.description=SEQ.id
		extractedRegions.append(SEQ)
		print(SEQ)
	
	with open(a.o+'regions.fasta', "w") as output_handle:
		SeqIO.write(extractedRegions, output_handle, 'fasta')

	if len(extractedRegions)>0:
		## Run ITSx on regions.fasta
		terminal('ITSx -t F -i %s -o %s --cpu %s --save_regions %s --not_found F --graphical F --summary F --positions F --fasta F' % (a.o+'regions.fasta',a.o+name,a.cpu,a.which))

		## Remove Duplicates from ITSx output
		with open(a.o+name+'.'+a.which+'.fasta','r') as handle:
			its1s=list(SeqIO.parse(handle,'fasta'))
		seqs=[str(s.seq) for s in its1s]
		diff=list(set(seqs))
		unique=[]
		for n in range(len(diff)):
			if len(diff)==1:
				Name=name
			else:
				Name=name+'_#'+str(n+1)
			unique.append(SeqRecord(Seq(diff[n]),Name,Name,str(seqs.count(diff[n]))+' copies of this sequence found in '+a.i.split('/')[-1]))
	
		with open(a.o+name+'.'+a.which+'_filtered.fasta', "w") as output_handle:
			SeqIO.write(unique, output_handle, 'fasta')

		## End message

		if len(diff)==1:
			print('\nDONE. A single '+a.which+' sequence was found in '+str(seqs.count(diff[n]))+' copies, in file '+a.i.split('/')[-1]+'.')
		elif len(diff)>1:
			print('\nDONE. '+str(len(diff))+' different '+a.which+' sequences were found in file '+a.i.split('/')[-1]+'. See output files for more info.')
		elif len(diff)==0:
			print('\nDONE. No '+a.which+' sequence found in file '+a.i.split('/')[-1])

	else:
		print('\nDONE. No '+a.which+' sequence found in file '+a.i.split('/')[-1]+'.')
		terminal('touch '+a.o+name+'.'+a.which+'_filtered.fasta')
		

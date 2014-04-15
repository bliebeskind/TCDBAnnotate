#! /usr/bin/env python

## For cleaning fasta description lines while retaining unique ids.
##
## Usage: clean_fasta_lines.py <infile> <outfile>

import sys
from Bio import SeqIO

def clean(infile,format='fasta'):
	'''
	Clean fasta description lines and return iterator over cleaned recs.
	
	Retains everything left of the first space unless the records is from JGI, 
	in which case the rec is is the second field (organism id) plus the second 
	field (protein id).
	'''
	records = SeqIO.parse(infile,format)
	for rec in records:
		if rec.id.startswith("jgi"):
			split = rec.id.split("|")
			rec.id = split[1] + split[2]
			rec.description = ''
		else:
			rec.id = rec.id.split(" ")[0]
			rec.description = ''
		yield rec
		
def write(infile,outfile,format='fasta'):
	'''Given infile, write outfile with cleaned description lines.'''
	cleaned_recs = clean(infile,format)
	SeqIO.write(cleaned_recs,outfile,'fasta')
	
if __name__ == '__main__':
	Usage = '\nUsage: clean_fasta_lines.py <infile> <outfile>\n'
	try:
		infile = sys.argv[1]
		outfile = sys.argv[2]
		write(infile,outfile)
	except IndexError:
		print Usage

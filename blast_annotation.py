#! /usr/bin/env python

import sys, pickle
import pandas as pd
from Bio.Blast import NCBIXML

## To do: 	Create functions for different database description line formats.
##			Integrate with sql3

## Ideas about sqlite schema:
## 	Fasta table: transcript name, nuc sequence, prot sequence, length
## 	TCDB table: transcript name, uniprot id, tcdb family id, tcdb family name, evalue
##	Uniprto table:
## 	HMMER table: 

## Also, might be good to have make_mapping be a generator, and have some other
## function to collect into memory and dump to pickle that is not used by default.

# Titles when blasting the TCDB look like this:
# gnl|BL_ORD_ID|125 gnl|TC-DB|Q9HCF6 1.A.4.5.6 Transient receptor potential cation channel subfamily M member 3 OS=Homo sapiens GN=TRPM3 PE=2 SV=4

def hit_generator(INFILE, EVALUE_THRESH):
	'''Generator function yielding tuple of query name and name of top hit.'''
	with open(INFILE) as f:
		blast_records = NCBIXML.parse(f)
		for record in blast_records:
			try:
				tophit = record.descriptions[0]
			except IndexError: # no hits?
				continue
			if float(tophit.e) < float(EVALUE_THRESH):
				yield record.query, tophit.title
			else:
				continue
				
def make_mapping(infile, evalue):
	hit_counter = 0
	bad_titles = []
	D = {"Query":[], "Uniprot Id":[], "TCDB Id":[], "Title":[]}
	print "\t".join(["Query","Uniprot Id","TCDB Id","Title"]) #header
	for query, title in hit_generator(infile, evalue):
		hit_counter +=1
		if hit_counter % 100 == 0:
			sys.stderr.write(str(hit_counter)+'\n')
		try:
			info = title.split(" ")
			ids = info[1].split('|')
			uniprot, tcdb, title = ids[2],ids[3],' '.join(info[2:])
		except (IndexError,AssertionError):
			bad_titles.append(title)
			continue
		D["Query"].append(query)
		D["Uniprot Id"].append(uniprot)
		D["TCDB Id"].append(tcdb)
		D["Title"].append(title)
		print "\t".join([query, uniprot, tcdb, title])
	sys.stderr.write("Found %i sequences with bad titles\n" % len(bad_titles))
	return D,bad_titles
		
			
if __name__ == '__main__':
	INFILE = sys.argv[1]
	EVALUE_THRESH = float(sys.argv[2])
	pickle_file = sys.argv[3]
	D,bad_titles = make_mapping(INFILE, EVALUE_THRESH)
	df = pd.DataFrame(D)
	with open(pickle_file,'w') as f:
		pickle.dump(df,f)
	with open("bad_titles.txt",'a') as f:
		for title in bad_titles:
			f.write(title)

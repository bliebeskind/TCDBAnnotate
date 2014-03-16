#! /usr/bin/env python

import sys, re
from pprint import pprint as pp

def parse(infile, cutoff):
	'''Return dictionary of hits better than cutoff mapped to e-values and 
	number of predicted domains, as a tuple.'''
	D = {}
	with open(infile) as f:
		for line in f:
			if not line.startswith('#'):
				line = line.strip().split(' ')
				linelist = [i for i in line if i != '']
				target, evalue, ndom = linelist[0],linelist[4], linelist[11]
				if float(evalue) < cutoff:
					D[target] = (evalue, ndom)
	return D

def file_list_from_infile(infile):
	'''Return list of file names from infile made with ls *.txt >> infile.
	Infiles must be of format species_* for use with function below.'''
	with open(infile) as f:
		return [line.strip() for line in f]
	
def hits_D_from_files(file_list, cutoff):
	'''Make a dictionary with species names as keys and dictionaries from
	parse() as values. Infiles must be of format species_*.'''
	D = {}
	for f in file_list:
		hits_D = parse(f, cutoff)
		species = f.split('_')[0] # Generate species name from infile name
		D[species] = hits_D
	return D

### Summary function.

def filtered_keys(key_list):
	expr = re.compile("T\d")
	nonredlist, redlist, genes = [], [], []
	try:
		for key in key_list:
			match = re.search(expr,key)
			assert match # did it find a match?
			assert match.end() == len(key) # was match at end of string?
			gene = key[:match.start()]
			if gene in genes:
				redlist.append(gene)
			else:
				nonredlist.append(key)
				genes.append(gene)
	except AssertionError: #not from Origins database
		return key_list, [] # return keylist as is and a blank redlist
	return nonredlist, redlist
	

def num_domains(domain_list):
	'''Print distribution of domains numbers.'''
	domain_set = set(domain_list)
	for item in domain_set:
		print "    Proteins with %s domain(s): %i"\
		% (item, domain_list.count(item))

def print_summary(D):
	'''Summarize hits with average e-value and distribution over number of 
	domains in the set of proteins..'''
	origins_redlist = []
	for key in D.iterkeys():
		species = D[key]
		print key+':'
		species_hits, redlist = filtered_keys(species.keys())
		origins_redlist += redlist
		numHits = len(species_hits)
		domain_list = [species[k][1] for k in species_hits]
		evalue_list = [float(species[k][0]) for k in species_hits]
		print "    Average e-value: %e" % (sum(evalue_list)/len(evalue_list))
		num_domains(domain_list)
		print ''
	print "** Origins of Multicellularity genes with multiple transcripts **"
	pp(origins_redlist)
	
		
if __name__ == '__main__':
	import cPickle as pickle
	infile, outfile = sys.argv[1], sys.argv[2]
	file_list = file_list_from_infile(infile)
	D = hits_D_from_files(file_list, .01)
	print_summary(D)
	pickle.dump(D, open(outfile, "wb"))

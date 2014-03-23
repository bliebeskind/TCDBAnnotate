#! /usr/bin/env python

import pandas as pd
import sys
from Bio import SeqIO


def get_genes(infile, gene_family):
	'''
	Infile is tab separated tsv files made by annotate.py,
	and then separated by TCDBAnnotate.to_tsv. Return dictionary of gene 
	names (first column) mapped to the name of their gene family.
	'''
	df = pd.read_table(infile,sep='\t')
	return {i: str(gene_family) for i in df["Query"]}
	
def make_gene_D(gene_list):
	'''Given list of gene family names, call gene_D and concatenate
	dictionaries.'''
	gene_D = {}
	family_count = 0
	for d in gene_list:
		try:
			with open(d+'.tsv') as f:
				gene_D.update(get_genes(f,d))
				family_count +=1
		except IOError:
			sys.stderr.write("No file called %s\n" % d)	
	print "Found %i gene families" % family_count
	return gene_D
	
def make_out_D(infile, gene_list):
	'''Given a fasta infile, and a list of gene family names, make
	dictionary with family names as keys and lists of SeqRecords
	as values. Call make_gene_D to get schema separating records
	into families.'''
	records = SeqIO.parse(infile,'fasta')
	gene_D = make_gene_D(gene_list)
	out_D = {}
	found = 0
	not_found = 0
	for rec in records:
		try:
			family = gene_D[str(rec.id)]
			if family not in out_D:
				out_D[family] = [rec]
			else:
				out_D[family].append(rec)
			found +=1
		except KeyError:
			not_found +=1
	sys.stderr.write("Found %i genes, didn't find %i\n" % (found, not_found))
	return out_D


def make_fasta_files(infile, gene_list):
	'''Take fasta infile and separate sequences into different families
	in gene list. Must have separate tsv files in directory made by 
	TCDBAnnotate.to_tsv'''
	out_D = make_out_D(infile, gene_list)
	total_seqs = 0
	total_families = 0
	for family in out_D:
		seq_list = out_D[family]
		num_seqs = SeqIO.write(seq_list,family+'.fas','fasta')
		total_seqs += num_seqs
		total_families += 1
	sys.stderr.write("Wrote %i sequences from %i families\n" %\
	(total_seqs, total_families))
	

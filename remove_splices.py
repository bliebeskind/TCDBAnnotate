#! /usr/bin/env python

from Bio import SeqIO
import sys

def get_longest(gene_D):
	'''Given a dictionary mapping gene names to dictionaries of their protein 
	products, each mapped to their length, return longest gene name'''
	func = lambda x,y,D: x if D[x] >= D[y] else y
	for gene,prots in gene_D.iteritems():
		func = lambda x,y: x if prots[x] >= prots[y] else y
		longest = reduce(func,prots.keys())
		return longest
	
def longest_variants_ensembl(infile,format='fasta'):
	'''Return generator of longest proteins from an Ensembl proteome. Will
	hit an AssertionError if all transcripts from each gene are not together
	in the proteome.'''
	proteins = SeqIO.parse(infile,format)
	all_genes = []
	gene_map = {} 	# {gene1: {prot1:24, prot2:109}}
	record_map = {} # {seq.description: SeqRecord}
	for rec in proteins:
		gene = rec.description.split(' ')[3]
		prot_len = len(str(rec.seq))
		if gene in gene_map:
			assert rec.description not in gene_map[gene],\
			"Protein name %s found twice" % rec.description
			gene_map[gene][rec.description] = prot_len
			record_map[rec.description] = rec
		elif gene_map: # gene map is populated
			assert gene not in all_genes, "Genes not continuous: %s" % gene
			all_genes.append(gene)
			longest = get_longest(gene_map) # get longest prot from previous 
			yield record_map[longest]		# gene and yield it
			gene_map,record_map = {},{} # flush dicts
			gene_map[gene] = {rec.description:prot_len}
			record_map[rec.description] = rec
		else: # gene map is not populated, i.e. first iteration
			gene_map[gene] = {rec.description:prot_len}
			record_map[rec.description] = rec
			all_genes.append(gene)
	longest = get_longest(gene_map) # get longest prot from last gene
	yield record_map[longest]		# and yield it
	sys.stderr.write("Found %i genes\n" % len(all_genes))
	
def longest_variants_origins(infile,format='fasta'):
	'''Return generator of longest proteins from an Origins of Multicellularity
	proteome. Will hit an AssertionError if all transcripts from each gene are 
	not together in the proteome.'''
	proteins = SeqIO.parse(infile,format)
	all_genes = []
	gene_map = {} 	# {gene1: {prot1:24, prot2:109}}
	record_map = {} # {seq.description: SeqRecord}
	for rec in proteins:
		gene = rec.description.split('|')[1].strip()
		prot_len = len(str(rec.seq))
		if gene in gene_map:
			assert rec.description not in gene_map[gene],\
			"Protein name %s found twice" % rec.description
			gene_map[gene][rec.description] = prot_len
			record_map[rec.description] = rec
		elif gene_map: # gene map is populated
			assert gene not in all_genes, "Genes not continuous: %s" % gene
			all_genes.append(gene)
			longest = get_longest(gene_map) # get longest prot from previous 
			yield record_map[longest]		# gene and yield it
			gene_map,record_map = {},{} # flush dicts
			gene_map[gene] = {rec.description:prot_len}
			record_map[rec.description] = rec
		else: # gene map is not populated, i.e. first iteration
			gene_map[gene] = {rec.description:prot_len}
			record_map[rec.description] = rec
			all_genes.append(gene)
	longest = get_longest(gene_map) # get longest prot from last gene
	yield record_map[longest]		# and yield it
	sys.stderr.write("Found %i genes\n" % len(all_genes))
		
if __name__ == '__main__':
	infile = sys.argv[1]
	outfile = sys.argv[2]
	database = sys.argv[3]
	if database == 'ensembl':
		SeqIO.write(longest_variants_ensembl(infile,format='fasta'),outfile,'fasta')
	elif database == 'origins':
		SeqIO.write(longest_variants_origins(infile,format='fasta'),outfile,'fasta')
	else:
		raise Exception("Database not found")

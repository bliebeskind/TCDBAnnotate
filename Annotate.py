#! /usr/bin/env python

import sys, os
import fasta_channel_families
import pandas as pd
import csv

## To do: 
## Would be nice not to have to write tsvs before fastas (see notes below)
## Could modify to_tsv to take in pickle files optionally, since 
## 	blast_annotation writes these.

tcdb_id_map = {"1.A.1.1": "KcsA",
		"1.A.1.2": "Kv",
		"1.A.1.3": "Slo",
		"1.A.1.4": "Kplant",
		"1.A.1.5": "CNG-HCN",
		"1.A.1.6": "MthK",
		"1.A.1.7": "Ktp-plant",
		"1.A.1.8":  "Ktp-animal",
		"1.A.1.9": "Ktp-animal",
		"1.A.1.10": "Nav",
		"1.A.1.11.1": "Cav",
		"1.A.1.11.2": "Cav",
		"1.A.1.11.3": "Cav",
		"1.A.1.11.4": "Cav",
		"1.A.1.11.5": "Cav",
		"1.A.1.11.6": "Cav",
		"1.A.1.11.7": "Cav",
		"1.A.1.11.8": "Cav",
		"1.A.1.11.9": "Cav",
		"1.A.1.11.10": "Leak",
		"1.A.1.11.11": "Cav",
		"1.A.1.11.12": "Cav",
		"1.A.1.11.13": "TPC",
		"1.A.1.11.14": "Cav",
		"1.A.1.11.15": "Leak",
		"1.A.1.11.16": "Leak",
		"1.A.1.11.17": "Leak",
		"1.A.1.11.18": "TPC",
		"1.A.1.11.19": "TPC",
		"1.A.1.11.20": "Cav",
		"1.A.1.11.21": "Cav",
		"1.A.1.11.22": "TPC",
		"1.A.1.11.23": "Leak",
		"1.A.1.11.24": "Cav",
		"1.A.1.12": "Kvirus",
		"1.A.1.13": "Kbac",
		"1.A.1.14": "BacNav",
		"1.A.1.15": "KCNQ",
		"1.A.1.16": "Sk",
		"1.A.1.17": "KvAP",
		"1.A.1.18": "TRESK",
		"1.A.1.19": "CatSper",
		"1.A.1.20": "Erg",
		"1.A.1.21": "NaK",
		"1.A.1.22": "MmaK",
		"1.A.1.23": "Kroot",
		"1.A.1.24": "CNR",
		"1.A.1.25": "MlotiK",
		"1.A.1.26": "Kplasmodium",
		"1.A.1.27": "Kputative1",
		"1.A.1.28": "Kputative2",
		"1.A.1.29": "Kputative3",
		}

class TCDBAnnotate:
	'''
	For handling sequences annotated by blast_annotation.py using the TCDB
	as a database. Can separate annotation file by gene and write out
	separate tsv and fasta files.
	
	Examples:
		>>> from Annotate import TCDBAnnotate
		>>> a = TCDBAnnotate()
	
	  ## Look at available gene names ##
		>>> a.family names()
	
	  ## Genes to look for are held in field "genes" ##
	  	>>> a.genes
	  	[]
	  	>>> a.append_genes(a.family_names()) # append all possible gene families
	  	>>> a.append_genes(["Kv","Nav'] # ...or only a few
	
	  ## Write out tsv and fasta files for all appended genes ##
		>>> a.to_tsv("my_annotation.tsv") # infile is blast_annotation output
		>>> a.to_fasta("my_annotation.tsv","my_fasta.fas")
		>>> ## Where 'my_fasta' is query file used for annotation.
	'''
	
	def __init__(self,schema=tcdb_id_map):
		self.schema = schema
		self.schema_genes = sorted(list(set([gene for gene in self.schema.itervalues()])))
		self.genes = []
		
	def family_names(self):
		'''Returns a list of gene names in schema.'''
		return self.schema_genes
		
	def append_genes(self,gene_list):
		'''Append gene_list to field value "genes". Each call will clear
		whatever is currently in this field.'''
		self.genes = [] # clear genes
		new_genes = []
		for gene in gene_list:
			try:
				assert gene in self.schema_genes
				new_genes.append(gene)
			except AssertionError:
				print "%s not in schema names" % (gene)
		self.genes += new_genes
	
	#not currently being used by anything
	def _make_directory(self,key):
		'''Make a directory named 'key' to write into'''
		if not os.path.exists("./%s" % key):
			os.mkdir("./%s" % key)
		else:
			sys.stderr.write("Directory %s already exists" % key)
			raise
	
	def _is_cav(self,tcdb_id):
		'''Check if tcdb_id belongs to the Cav group, which unfortunately
		includes NALCN, TPCs, and others.'''
		if "1.A.1.11." in tcdb_id:
			return True
		else:
			return False
		
	def _gene_generator(self,infile,gene_list):
		'''Generator function for parsing tsv annotation file.
		Infile: *tsv* output from blast_annotation.py.'''
		with open(infile) as f:
			for line in f:
				line = line.split("\t") # must be tab-separated
				assert len(line) == 4, "Bad annotation file - must be tsv."
				if line[0] == "Query": # skip header
					continue
				tcdb_id = line[2]
				if self._is_cav(tcdb_id):
					pass
				else:
					tcdb_id = ".".join(tcdb_id.split(".")[:4])
				try:
					if self.schema[tcdb_id] in gene_list:
						yield "\t".join(line), self.schema[tcdb_id]
				except KeyError:
					pass # for now...
					# schema should eventually include everything
					# and use this code:
					#raise Exception("%s not found in schema" % tcdb_id)
					
	def _gene_annotation_D(self,infile,gene_list):
		'''Call _gene_generator and make a dictionary that maps gene names
		to list of lines from original annotation file that have these genes.'''
		gene_D = {gene:["\t".join(["Query","Uniprot Id","TCDB Id","Title\n"])]\
			for gene in gene_list}
		for line,gene in self._gene_generator(infile,gene_list):
			assert gene in gene_D, "Parsing problem: %s" % gene
			gene_D[gene].append(line)
		return gene_D
					
	def to_tsv(self, infile, gene_list=None):
		'''Write separate tsv files for each gene in gene_list. 
		Infile: a *tsv* file output by blast_annotation.py
		'''
		if gene_list == None:
			gene_list = self.genes
		gene_D = self._gene_annotation_D(infile,gene_list)
		for gene,L in gene_D.iteritems():
			if len(L) == 1: # no hits found for gene (only header line)
				continue
			with open(gene+'.tsv','w') as outfile:
				outfile.write("".join(L))
			
	# Currently works but would be nice not to have to write
	# tsv files first. Could maybe pickle the information?
	def to_fasta(self,fasta_file,gene_list=None):
		'''Take a fasta file that of queries that were annotated and separate
		into fasta files by gene in gene_list.'''
		if gene_list == None:
			gene_list = self.genes
		fasta_channel_families.make_fasta_files(fasta_file,gene_list)


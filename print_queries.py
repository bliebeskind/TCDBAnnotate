#! /usr/bin/env python

## This can probably be expanded to be more flexible and be able to 
# print whatever info is put into the script.

import blast_parser
import sys

def blast_info_list(searches):
    '''Returns list of blast info dictionaries from a parsed blast hits
    file.'''
    blast_info = []
    for search in searches:
        blast_info.append(search[0])
    return blast_info
        
def get_queries(searches):
    '''Returns list of queries from list of blast info dictionaries.'''
    search_counter = 0
    queries = []
    for info in blast_info_list(searches):
        search_counter += 1
        try:
            queries.append(info['Query'].strip())
        except KeyError:
            sys.stderr.write('Search %i has no query' % search_counter)
    return queries, search_counter

def print_queries(filename):
    '''Prints all queries from a Blast hits file.'''
    searches = blast_parser.get_searches_from_file(filename)
    queries, search_counter = get_queries(searches)
    for query in queries:
        print query   #prints queries even if redundant - fix?
    print '\nPrinted queries from %i searches' % search_counter
    
if __name__ == '__main__':
    Usage = "print_queries <file>"
    try:
        print_queries(sys.argv[1])
    except IndexError:
        print Usage

#!/usr/bin/python2.6
#
# filters a fasta file (from standard input), pulling out sequences by either taxid or taxname, using
# NCBI tax files - see below. Use-case is to pull out seqs from a specific (older) fasta and 
# taxid file combination  
#
# Requires biopython
#
# Alan McCulloch 2/2014
#
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import sys
import os
import itertools
import re
import string
 
usage = """
usage :  cat some_file.fa | tax_yank.py -byid|-byname tax_dump_file [taxnamefile] tax1,tax2,tax3 ....

tax_dump_file is an NCBI tax dump file containing 2 columns
gi taxid

tax1 tax2 etc are either comma-seperated taxid numbers, or comma seperate species names 

taxnamefile is an NCBI tax file giving names and descriptions of all taxid numbers

-e.g. 

1       |       all     |               |       synonym |
1       |       root    |               |       scientific name |
2       |       Bacteria        |       Bacteria <prokaryote>   |       scientific name |
2       |       Monera  |       Monera <Bacteria>       |       in-part |
2       |       Procaryotae     |       Procaryotae <Bacteria>  |       in-part |


The script will filter the sequence data, only outputting sequences with gi numbers 
that are associated with one of tax1 tax2 etc (via the mappings from tax_dump_file (and optionally tax name file))

Examples : 

from /localdata/mascotdb2012/taxonomy/names.dmp

9031    |       Gallus gallus   |               |       scientific name |
562     |       Escherichia coli        |               |       scientific name |
9913    |       Bos taurus      |               |       scientific name |

# this gives a smallish file , only including seqs with three specific taxid
cat /localdata/mascotdb2012/sequence/NCBInr/current/NCBInr_20121217.fasta |  ./tax_yank.py -byid /localdata/mascotdb2012/taxonomy/gi_taxid_prot.dmp 9031,562,9913

# this gives a larger file, since we use all taxid that match the given names
cat /localdata/mascotdb2012/sequence/NCBInr/current/NCBInr_20121217.fasta |  ./tax_yank.py -byname /localdata/mascotdb2012/taxonomy/gi_taxid_prot.dmp /localdata/mascotdb2012/taxonomy/names.dmp Gallus gallus,Escherichia coli,Bos taurus

"""

TAXOPTION_ARGINDEX=1
TAXFILE_ARGINDEX=2
TAXIDS_ARGINDEX=3
TAXNAMEFILE_ARGINDEX=None

TAXID_INDEX=1
SEQGI_INDEX=1

 
if len(sys.argv) <  TAXIDS_ARGINDEX:
   print usage
   sys.exit(1)

if sys.argv[TAXOPTION_ARGINDEX] not in ["-byid", "-byname"]:
   print usage
   sys.exit(1)

if sys.argv[TAXOPTION_ARGINDEX] == "-byname":
   TAXNAMEFILE_ARGINDEX = 3
   TAXIDS_ARGINDEX = 4

if not os.path.isfile(sys.argv[TAXFILE_ARGINDEX]):
   print "%s is not a file"%(sys.argv[TAXFILE_ARGINDEX])

if TAXNAMEFILE_ARGINDEX is not None:
   if not os.path.isfile(sys.argv[TAXNAMEFILE_ARGINDEX]):
      print "%s is not a file"%(sys.argv[TAXNAMEFILE_ARGINDEX])

# if necessary get numeric taxids from tax names
taxidlist=re.split("\s*\,\s*", string.join(sys.argv[TAXIDS_ARGINDEX:], " "))
taxiddict=dict(zip(taxidlist, taxidlist))
if sys.argv[TAXOPTION_ARGINDEX] == "-byname":
   # get the names passed in as args, into a regular expression like (species1)|(species2)| etc
   nameRegexp = string.join(["(%s)"%species for species in taxidlist],"|")

   # get an iterator through the names file
   taxnameiter = (re.split("\s*\|\s*", record.strip())[0:2] for record in open(sys.argv[TAXNAMEFILE_ARGINDEX], "r"))
   # iterates records like ['9', 'Buchnera aphidicola Munson et al. 1991']
      
   # make a filtered iterator for records that match the regexp
   selectednameiter = itertools.ifilter(lambda x:re.search(nameRegexp, x[1], re.IGNORECASE) is not None, taxnameiter)
   taxiddict = dict(selectednameiter) # a dictionary with key = taxid

# now get gilist
# raw iterator : 
taxdumpiter = (re.split("\s+",record.strip()) for record in open(sys.argv[TAXFILE_ARGINDEX],"r"))
# filtered iterator:
taxiter = itertools.ifilter(lambda x:x[TAXID_INDEX] in taxiddict , taxdumpiter)
# make a gi lookup dictionary 
#taxdict = dict(itertools.islice(taxiter,0,5)) # for a shorter test
taxdict = dict(taxiter)

 
# each record assumed to have an id like gi|67475154|ref|XP_653293.1| - write to stdout records where gi is in our list
SeqIO.write ((r for r in SeqIO.parse(sys.stdin, "fasta") if re.split("\|", r.id)[SEQGI_INDEX] in taxdict) , sys.stdout, "fasta")
 

#!/usr/bin/env python

#########################################################################
#This program is free software: you can redistribute it and/or modify   #
#it under the terms of the GNU General Public License as published by   #
#the Free Software Foundation, either version 3 of the License, or      #
#(at your option) any later version.                                    #
#                                                                       #
#This program is distributed in the hope that it will be useful,        #
#but WITHOUT ANY WARRANTY; without even the implied warranty of         #
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
#GNU General Public License for more details.                           #
#                                                                       #
#                                                                       #
#You should have received a copy of the GNU General Public License      #
#along with this program.  If not, see <http://www.gnu.org/licenses/>.  #
#                                                                       #
# contact : alan.mcculloch@agresearch.co.nz                             #
#########################################################################
from types import *
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import ProtParam
import string
import itertools
import collections
import re
import sys
import exceptions
import os
import argparse


class myException(exceptions.Exception):
    def __init__(self,args=None):
        super(myException, self).__init__(args)


os.linesep = '\n' # for a Unix file
TBlastHit = collections.namedtuple('TBlastHit', ['qseqid' , 'qstart' , 'qend', 'qframe' , 'sacc', 'stitle'])
THitGroup  = collections.namedtuple('THitGroup', ['qseqid' , 'hitlist'])
TOrf = collections.namedtuple('TOrf', ['frame' , 'dna_start' , 'dna_end', 'seqrecord'])

def main():

    args  = getParameters()
 
    print "using %s\n\n\n"%str(args)

    getORFS(infile=args.infile,outfile=args.outfile, markup=args.markup, minlength = args.minlength, hitfile = args.blastfile, filter=args.filter)


def getParameters():
    description = """
          This script translates DNA to protein, optionally using blast hits to
          determine the translation frame """
    long_description = """
          This script translates DNA to protein, optionally using blast hits to
          determine the translation frame - if all hits are in the same frame ,
          then that frame is used. If there are hits in different frames, blast
          evidence is not used and the sequence is translated in all 6 frames.
          If there are no hits, the sequence is translated in all 6 frames).

          It is designed for processing large files, so uses iterators to consume
          streams of DNA sequences and hits, which are assumed to be in the same
          order (e.g. rather than loading a file of hits into a dictionary).
          Hits can be piped to the script (for example via gunzip -c myhits.gz, or
          a blast command can be piped directly to the script, avoiding the need
          for storing a file of hits)

          If an output filename is not specified , translations are written to stdout

          If no blast hits are supplied, all translations are in all frames

          Examples :

          # use a file of hits
          big_blast_translate.py Transcript_and_mRNA-seq_contigs.fa -o Transcript.translated -b Transcripts.hits

          # pipe hits from a compressed file to the program
          gunzip -c Transcripts.hits.gz   |  big_blast_translate.py Transcript_and_mRNA-seq_contigs.fa -o Transcript.translated -b -

          # pipe blast results to the program, avoiding the need to store a file of hits
          blastx -query Transcript_and_mRNA-seq_contigs.fa -num_threads 4 -db uniprot_sprot.fa -evalue 1.0e-6 -outfmt "7 qseqid qstart qend qframe stitle" | big_blast_translate.py Transcript_and_mRNA-seq_contigs.fa -o Transcript.translated -b -

          # piping blast to the this script as above but run in chunks on a cluster using tardis
          tardis.py -w -c 1000 -d my_scratch_dir blastx -query _condition_fasta_input_Transcript_and_mRNA-seq_contigs.fa -num_threads 4 -db uniprot_sprot.fa -evalue 1.0e-6 -outfmt \"7 qseqid qstart qend qframe stitle\" \| big_blast_translate.py -i _condition_fasta_input_Transcript_and_mRNA-seq_contigs.fa -o _condition_text_output_Transcript.translated -b -

          # don't use any hits - all seqs translated in 6 frames
          big_blast_translate.py -i Transcript_and_mRNA-seq_contigs.fa -o Transcript.translated

          # filter the output - use a regular expression which is applied to the description text of each blast hit
          big_blast_translate.py -i Transcript_and_mRNA-seq_contigs.fa -o Transcript.translated -f "GN\=\S+"

          The hitfile is assumed to be formatted as per the output of a blast run like the following :

          blastx -query inseqs.fa -num_threads 4 -db uniprot_sprot.fa -evalue 1.0e-6 -outfmt "7 qseqid qstart qend qframe sacc stitle"

          which contains records like

# BLASTX 2.2.28+
# Query: 140304CS1900602900001 CS1900602900001 DNA=(0%N 18%C 19%G 31%A 32%T ) Expr= 4=CANT(1) CEHF(1) CEFE(1) CEEF(1) (no refseq hit) (no nr protein hit)
# Database: uniprot_sprot.fa
# 0 hits found
# BLASTX 2.2.28+
# Query: 140304CS19010313FFFFB CS19010313FFFFB DNA=(0%N 21%C 25%G 29%A 24%T ) Expr= 1=CEIE(1) Refseq=gi|77404355|ref|NM_005801.3| Homo sapiens eukaryotic translation initiation factor 1 (EIF1), mRN(eval= .00
E+00) NRProtein=ref|XP_418159.2| PREDICTED: similar to isolog of yeast sui1 and rice gos2; putative [Gallus gallus(eval= 2.00E-45)
# Database: uniprot_sprot.fa
# Fields: query id, q. start, q. end, query frame, subject acc., subject title
# 23 hits found
140304CS19010313FFFFB   189     500     3       Q5RFF4  Eukaryotic translation initiation factor 1 OS=Pongo abelii GN=EIF1 PE=3 SV=1
140304CS19010313FFFFB   189     500     3       P41567  Eukaryotic translation initiation factor 1 OS=Homo sapiens GN=EIF1 PE=1 SV=1
140304CS19010313FFFFB   189     500     3       Q5E938  Eukaryotic translation initiation factor 1 OS=Bos taurus GN=EIF1 PE=3 SV=1
140304CS19010313FFFFB   189     500     3       P48024  Eukaryotic translation initiation factor 1 OS=Mus musculus GN=Eif1 PE=2 SV=2
140304CS19010313FFFFB   189     500     3       P61220  Eukaryotic translation initiation factor 1b OS=Sus scrofa GN=EIF1B PE=3 SV=1
140304CS19010313FFFFB   189     500     3       Q9CXU9  Eukaryotic translation initiation factor 1b OS=Mus musculus GN=Eif1b PE=2 SV=2
140304CS19010313FFFFB   189     500     3       Q4R4X9  Eukaryotic translation initiation factor 1b OS=Macaca fascicularis GN=EIF1B PE=3 SV=1
140304CS19010313FFFFB   189     500     3       O60739  Eukaryotic translation initiation factor 1b OS=Homo sapiens GN=EIF1B PE=1 SV=2
etc

    """
    parser = argparse.ArgumentParser(description=description, epilog=long_description, formatter_class = argparse.RawDescriptionHelpFormatter)
    parser.add_argument("infile", help="name of input fasta file of DNA to translate")
    parser.add_argument("-o","--outfile", help="name of output file (default stdout)")
    parser.add_argument("-b","--blastfile", help="name of file containing blast hits used to determine frame (- = stdin)")
    parser.add_argument("-f","--filter", help="regular expression to filter sequences, applied to sequence description")
    parser.add_argument("-l","--minlength", help="minimum length of translated sequence to output (default 20)" , type=int, default=20)
    parser.add_argument("-m","--markup", help="markup descriptor controlling annoation markup added to sequence description (default: frame,start,end,protparam)" , default="frame,start,end,protparam")

    args = parser.parse_args()
    return args

def merge(seqs, raw_hits):
    """
    does a "grouped left outer join" of seqs and hits - i.e. yields (seq [hit list]|None)
    for all seqs (assuming all hits have seqs, some seqs have hits, seqs have in general
    many hits, and hits are ordered in the file as are seqs)

    to use , call as (e.g.)

    for seqhit in merge(seqs, hits):
        print seqhit

    where seqs iterates over a fasta file
    raw_hits is a file stream, in which records are assumed to be provided in the format of tabular blast as above

    """

    # create an iterator that ignores the comment lines and parses hits into named tuples
    hits = itertools.ifilter(lambda record:record[0] != "#", raw_hits )
    hits = ( TBlastHit(recordtuple[0], int(recordtuple[1]), int(recordtuple[2]), int(recordtuple[3]), recordtuple[4], recordtuple[5]) for recordtuple in (re.split("\t", record.strip()) for record in hits if len(re.split("\t",record.strip())) == 6) )

    # group hits so for each sequence, we pull out a list of all its hits
    grouped_hits = itertools.groupby(hits, lambda x:x.qseqid)
    grouped_hits = iter((THitGroup(gh[0], list(gh[1])) for gh in grouped_hits))

    # implement the "left outer join" - i.e. for seqs with no hits , return
    # (seq, None)
    seq = seqs.next()

    have_hits = True

    try:
        hit_list = grouped_hits.next()
    except StopIteration:
        print >> sys.stderr, "*** Warning - no blast hits supplied, all translations will use all 6 frames ***"
        have_hits = False

    while have_hits:
        if seq.name == hit_list.qseqid:
            yield(seq,hit_list)
            seq = seqs.next()
            try:
                hit_list = grouped_hits.next()
            except StopIteration:
                have_hits = False
        else:
            yield(seq,None)
            seq = seqs.next()

    # after we have run out of hits, yield the rest of the seqs
    while True:
        yield(seq,None)
        seq = seqs.next()

def getORFS(infile, outfile,  markup , minlength, hitfile, filter):
    """
    Process an input file (FASTA format) and an optional input file of blast hits used to 
    determine consensus frame, write translations in either consensus frame or all 6 frames
    """
    fastaReader = SeqIO.parse(file(infile,"r") , "fasta") 

    if outfile != None:
        fastaWriter = file(outfile,"w")
    else:
        fastaWriter = sys.stdout

    # if hits have been supplied, get an iterator over seqs and hits  
    if hitfile == '-':
        seqs_with_hits = merge(fastaReader, sys.stdin)
    elif hitfile is not None:
        seqs_with_hits = merge(fastaReader, file(hitfile,"r"))
    else:
        seqs_with_hits = merge(fastaReader, iter(()))
 

    # for each sequence and any hits it has, decide a translation method, and call dnaToORFs to get translation(s)
                 
    recordCount = 0
    for (dnaRecord, hitgroup) in seqs_with_hits:
        recordCount += 1
        translationHint = ""
        if hitgroup is not None:
            hits = hitgroup.hitlist
            # count whether all hits have the same frame as the first

            if len(hits) != len([item for item in hits if item.qframe == hits[0].qframe]):
                translationHint="(translation method=6 frames because blast hits in different frames)"
                hit = None
            else:
                hit = hits[0] # assume hits are sorted in the file so the best is first
                translationHint="(translation method=using blast hit consensus frame %s)"%hit.qframe

        else:
            translationHint="(translation method=6 frames because no blast hits)"
            hit = None

        if filter is not None and hit is None:
            continue

        orfs = dnaToORFs(dnaRecord,minlength,markup, hit, translationHint, filter)
      
        SeqIO.write(orfs, fastaWriter, "fasta")

    fastaWriter.flush()
    fastaWriter.close()
    

def dnaToORFs(dnaRecord, minlength=20, markup=['frame','start','end'], hit = None, translationHint="", filter=None):
    """ this method returns a list of seq objects, each being an ORF from the translation
    of the given DNA.  
    """
    seq=dnaRecord.seq

    # lower case symbols are sometimes interpreted as stops - upper-case everything 
    seq = seq.upper()
 
    
    orfs=[]
    hitorfs = [] # used if we are to filter out orfs that match a given frame and overlap a hit
    frames = {
        1 : 'F1',
        2 : 'F2',
        3 : 'F3',
        -1 : 'R1',
        -2 : 'R2',
        -3 : 'R3'
    }

    # process translations in each frame
    for frame in frames.keys():
        # get translation with * for stops
        if frame > 0:
            translation = seq[frame-1:].translate()
        else:
            translation = seq.reverse_complement()[abs(frame)-1:].translate()
         
        # the translation gives, e.g. , DKYTHTHTHTHTKSLLCVRHCAKYWGYKSKE*NDPYS*GIYLLSGKTAYMYKHTE*IQSRERGC*E
        # - parse out ORFS using a split on  * character
        iter = translation.split("*")
        start =  abs(frame)  # this will be updated to equal start of each peptide in the orignal DNA   
                             # (for reverse strand frames, is the start in the reverse complement sequence)
        for ipeptide in range(0,len(iter)):
            peptide = iter[ipeptide]
            if len(peptide) < minlength:
                start = start + 3*len(peptide) + 3 # increment of 1 is due * seperator
                continue
            else:
                description = 'translation of %s %s '%(dnaRecord.name, translationHint)

                peptidestart = start # in case in future we start peptide at an offset - e.g. 
                                     # may want to add restriction that start AA must be M 
                # calculate end
                peptideend = peptidestart + 3*len(peptide) - 1

                # if frame < 0, then the position is relative to the start of the reverse-complement sequence. 
                # convert this so relative to start of orignal sequence
                if frame < 0:
                    peptidestart = 1 + len(seq) - peptidestart 
                    peptideend = 1 + len(seq) - peptideend


                # sanity check our start and stop against the length of the peptide
                if 1+abs(peptideend-peptidestart) != 3*len(peptide):
                    raise myException("Error for ORF in frame %s of %s , start=%s end=%s but peptide length = %s"%\
                                      (frame, dnaRecord.name, peptidestart, peptideend, len(peptide)))

                if 'frame' in markup:
                    description = "%s ; frame=%s"%(description, frame)

                if 'start' in markup:
                    description = "%s ; start=%s"%(description, peptidestart)

                if 'end' in markup:
                    description = "%s ; end=%s"%(description, peptideend)

                # do annotation of the peptide, if this has been requested
                analysis = None
                if 'protparam' in markup:
                    if not "X" in peptide:
                        analysis = ProtParam.ProteinAnalysis(str(peptide))
                        description = "%s ; molweight=%s"%(description, analysis.molecular_weight())
                        description = "%s ; aromaticity=%s"%(description, round(analysis.aromaticity(),4))
                        description = "%s ; instability_index=%s"%(description, round(analysis.instability_index(),4))
                        description = "%s ; gravy=%s"%(description, round(analysis.gravy(),4))
                        description = "%s ; isoelectric_point=%s"%(description, round(analysis.isoelectric_point(),4))
                        description = "%s ; secondary_structure_fraction=(%s)"%\
                                  (description, reduce(lambda x,y : x+','+y,[" %4.4f"%item for item in analysis.secondary_structure_fraction()]))
                if 'flexibility' in markup:
                    if analysis == None:
                        analysis = ProtParam.ProteinAnalysis(str(peptide))
                        
                    description = "%s ; flexibility=%s"%(description, analysis.flexibility())
                    
    
                # make a sequence record - this is the type expected by the BioPython  fasta writer class
                if hit != None:
                    description = hit.sacc + " " + hit.stitle + " : " + description
                    
                record = SeqRecord(\
                    peptide, \
                    id = "%s_%s_%s"%(dnaRecord.name,frames[frame],peptidestart), \
                    name = "%s_%s_%s"%(dnaRecord.name,frames[frame],peptidestart), \
                    description = description
                    )

                if hit != None:
                    if frame == hit.qframe:
                        if frame > 0:
                            hitorfs.append(   TOrf(frame, peptidestart, peptideend, record) )
                        else:
                            hitorfs.append(   TOrf(frame, peptideend, peptidestart, record) )

                            
                # append this to the list of ORFs we might return ( = will return, if we have no hit, also will if we have 
                # a hit but do not find any ORFS in the frame's hit (unusual))
                orfs.append(record)

                # accumulate the start , so will be correct for the next peptide in this scan
                start = start + 3*len(peptide) + 3 # increment of 3 is due * seperator

    # if the translation used a blast hit, pick the matching orf 
    if hit != None and len(hitorfs) > 0:
        # pick the longest ORF that matches the frame and overlaps the hit
        # if we don't find any overlaps - just return the longest one
        # first sort in descending order so that in case there are two we pick the longest
        hitorfs.sort(lambda x,y:abs((x.dna_end-x.dna_start)) - abs((y.dna_end-y.dna_start)) )
        hitorfs.reverse()

        best_orf = hitorfs[0].seqrecord

        for orf in hitorfs:
            # (note, if a negative frame, then qstart > qend , and dna_start > dna_end)
            if ( (hit.qstart <= orf.dna_start <= hit.qend)  or 
                 (hit.qend <= orf.dna_start <= hit.qstart) or
                 (orf.dna_start  <= hit.qstart <= orf.dna_end )  or
                 (orf.dna_end  <= hit.qstart <= orf.dna_start ) ): 

                best_orf = orf.seqrecord
                break

        if filter is not None:
            #print "searching for %s in %s"%(filter, best_orf.description)
            if re.search(filter, best_orf.description, re.IGNORECASE) is not None:
                orfs = [best_orf]
            else:
                orfs=[]
        else:
            orfs = [best_orf] 
    elif len(hitorfs) == 0 and filter is not None:
        orfs=[]

    #print "returning %s"%str(orfs)
    return orfs
            

if __name__ == "__main__" :
    main()

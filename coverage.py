#!/usr/bin/env python

import time
import resource
import sys
import os
from subprocess import Popen, PIPE
from optparse import OptionParser

parser = OptionParser()
parser.add_option("-1", "--1", dest="forward_reads",
                  help="forward reads file", metavar="FILE")
parser.add_option("-2", "--2", dest="reverse_reads",
                  help="reverse reads file", metavar="FILE")
parser.add_option("-x", "--x", dest="reference_fasta",
                  help="reference fasta", metavar="FASTA")
parser.add_option("-o", "--out", dest="outdir", metavar="STRING",
                  help="output dir", type="string", default="out")
parser.add_option("-q", "--q", dest="mean_quality_score", metavar="INT",
                  help="filter reads by mean quality score (range: 0-41)",
                  default=0, choices=[str(x) for x in range(0,41)])

def process_reads(quality, score_cutoff):
    reads_average_quality = 0
    for i in quality:
        reads_average_quality += ord(i)-33
    if float(reads_average_quality)/len(quality) > score_cutoff:
        return 1
    return 0

def filter_low_quality_alignment(infile, score_cutoff):
    readfile = open(infile, 'r')
    outfile = open(os.path.splitext(infile)[0] + '.filtered.sam', 'w')
    total_seq = 0
    qualified_seq = 0
    if_qualify = 0 # decision made based on filter criteria
    for line in readfile:
        if line[0] == '@':
            outfile.write(line)
            continue
        total_seq += 1
        items = line.rstrip().split('\t')
        seq = items[9]
        quality = items[10]
        if len(seq) != len(quality):
            raise Exception('Error in fastq: sequence length does no match quality score length')
        if process_reads(quality, int(score_cutoff)):
            outfile.write(line)
            qualified_seq +=1
    readfile.close()
    outfile.close()


if __name__ == "__main__":

    if sys.version < '2.7':
        raise Exception('Requires Python 2.7 or later.')

    start_time = time.time()

    soft, hard = 1.6*10**10, 1.6*10**10   ## 16G memory limit
    resource.setrlimit(resource.RLIMIT_AS, (soft, hard))

    (options, args) = parser.parse_args()
    if len(sys.argv) < 2:
        parser.print_help()
        exit(0)

    ##################
    ##     main     ##
    ##################

    try:
        ######### mkdir outfolder #########
        if not os.path.isdir(options.outdir):
            os.makedirs(options.outdir)

        working_dir = os.getcwd()

        ######### 1. alignment #########
        ## generate index file
        os.system('''module load bowtie2/2.2.3''')
        index_foler = os.path.dirname(os.path.abspath(options.reference_fasta))
        os.chdir(index_foler)
        os.system('''
            bowtie2-build {0} HIV
            '''.format(options.reference_fasta))
        os.chdir(working_dir)
        ## bowtie2 alignment
        os.system('''
            bowtie2 -p 1 -N 0 -x {0} -1 {1} -2 {2} -S {3}
            '''.format(index_foler+os.sep+'HIV', options.forward_reads, options.reverse_reads, options.outdir+os.sep+'align.sam'))

        ######### 2. filter low quality alignment #########
        filter_low_quality_alignment(options.outdir+os.sep+'align.sam', 
                                     options.mean_quality_score)

        ######### calculate coverage #########
        os.chdir(options.outdir)
        os.system('''module load samtools/1.4''')
        os.system('''samtools view -bS align.filtered.sam -o align.filtered.bam''')
        os.system('''samtools sort -o align.filtered.sort.bam align.filtered.bam''')
        os.system('''samtools depth -aa -d 1000000 align.filtered.sort.bam > align.filtered.sort.bam.depth''')
        ## generate coverage table and coverage plot
        os.system('''module load R/3.4.1-shlib''')
        os.system('''Rscript ../plot_coverage.R align.filtered.sort.bam.depth''')

        ######### 3. calculate correlation between coverage and seq content#########
        os.system('''samtools mpileup -Bf {0} -aa align.filtered.sort.bam > align.filtered.sort.bam.mpileup'''.format(options.reference_fasta))
        os.system('''Rscript ../calc_correlation.R align.filtered.sort.bam.depth align.filtered.sort.bam.mpileup''')
        os.system('''Rscript -e "rmarkdown::render('correlation.Rmd', 
								 params=list(nt_count={0}, nt_png={1}))"
				  '''.format(options.outdir+os.sep()+'coverage_nt.tsv', options.outdir+os.sep()+'avg_cov.png'))
    except MemoryError as err:
        sys.exit('memory exceeded')

    print("--- %s seconds ---" % (time.time() - start_time))
    

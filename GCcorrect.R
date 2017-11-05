
###
#
# Installing
#
###
# R CMD INSTALL GCcorrect_0.96.tar.gz


library("GCcorrect")

#####
# Begin by preparing a chromosome structure.
#
# reference: The reference chromosome,
#            Formats accepted:
#            1. A numeric vector with 4 values:
#               0,1,2,3,NA numeric (representing "A","T","C","G",NA)
#            2. A fasta file-name, or 0,1,2,3,NA numeric (representing "A","T","C","G",NA)
#
# repeats:   Unmappable positions in the chromosome
#            Formats accepted:
#            1. Logical vector in R format marking unmappable positions (0 - mappable, 1- repeat)
#            2. A fasta file-name similar to the reference, with "R" marking unmappable positions
#               (see example in GCcorrect/data/repeat75chr1.fa) 
#            We attach to files in GCcorrect/python/ to help create these files. 
#
#######

# The README file will give some background on this data...
args <- commandArgs(trailingOnly = TRUE)

chr1 = prepareChrom("HIV.RData",reference =args[1], repeats = args[2], details = list("HIV",1,200))

# If you prefer working with the full chromosome 1, you can find it in
# http://www.stat.berkeley.edu/share/yuvalb/HCCchr1.tar

# Reading the repeats file might give out a warning (file ended unexpectedly)
# That is fine. 

# Each chromosome should have a separate structure. 

#####
# Input reads
#
# A table for the forward strand and a table for reverse strand reads that were mapped to a chromosome. 
# Tables should have two columns:
# 1. 5' end of the read (as mapped to the forward strand)
# 2. Length of the fragment (if available, or 0 otherwise)
#    See comment below on length inputs; 
#    
# We usually keep chromosome as a first column, and cut it away. 
# See HCC1569_chr1_sm_for.cln and HCC1569_chr1_sm_rev.cln for the format 
#######

dat1_for = as.matrix(read.table(args[3])[,2:3])
dat1_rev = as.matrix(read.table(args[4])[,2:3])

####### Comments on Read Input    ##########
# Indexing starts at 1 (not 0).
# We assume index points to the 5' end of the read.
#
# Forward strand this will be the first base of the sequenced read.
# An easy way to check:
# chr1$chrline[dat1_for[1,1] + 0:9] should give first bases of the read.
# 
# Reverse strand index should point to the 5' of the complementary fragment.
# chr1$chrline[dat1_rev[1,1] + 0:9] should give first bases of the reverse complement read. 
# This is not(!) the end of the fragment.
# Some mapping software output this by default. Otherwise just subtract (read_len-1). 
#
# Note: Make sure not to use the same pair-end in both the forward and in the reverse strand.
#       An easy way to do this is to first run the matching of 5' and 3',
#       and then only use the 5' end results and the length information. 
#
# Length  
# We mean the number of bases between 5' end of fragment and the 3' end.
# This is not always the default of the mapping software
# (some exclude the second "read" for technical reasons)
#
# If length information is available, using only forward strand mappings
# would probably be sufficient. For single-stranded data though,
# splitting by strand is important.
#
################


#
####


# If you are using only forward strand, prepare a dummy variable for the reverse:
# dat1_rev = matrix(1,nc=2,nr=2);
# And make sure to always use strand="F" where relevant

dat1 = prepareReads (dat1_for, dat1_rev,chr_len = length(chr1$chrline), 1, 200,'HIV')


# Basic examination of the reads:
basicPlots(dat1,types = c(1,2,3,4))
# Plots types are:
# 1) Density of read counts
# 2) Length density curves - total, forward, reverse
# 3) Counts by genome location
# 4) Comparison of counts on between strands

# Examination of relation between reads and chromosome (mappability, GC)

compDataRefPlots(chr1 ,dat1,types = c(1,2,3,4))
# Plot types are:
# 1) Binned readcounts, mappability and GC across the chromosome
# 2) Plot of counts by mappability
# 4) Plot of GC by mappability
#
# 3) Calculates compares how likely a read is to be mapped in an unmappable location,
#    compared to the proportion of unmappable locations. This is another indication
#    of the quality of this mappability score

# To get reproduceable sample...
set.seed(1000)

# We want to create a manual mask for the genome,
# removing zero stretches and pileups
useSamp = !logical(length(chr1$isgc)) 
zeros_10K = findLongZeros(dat1$for10K,2)
tmp_size <- 50
if (!is.null(ncol(zeros_10K)) & ncol(zeros_10K) > 0) {
    for (i in 1:ncol(zeros_10K)) {
        useSamp[((tmp_size*(zeros_10K[1,i]-1)) + 1) : (tmp_size*(zeros_10K[2,i]))] = F
    }
}
useSamp = useSamp[1:length(chr1$isgc)]

high_10K = which(dat1$for10K>600)
for (i in high_10K) {
  useSamp[tmp_size*(i-1) + 1:tmp_size] = F
}
useSamp = useSamp[1:length(chr1$isgc)]

#####
# Create a sample of chromosome 1
#####
sampCh1 = sampleChrom(chr1,dat1,n=5000,margin = 1000,len_range = c(), memoryopt = TRUE,useSamp = useSamp) 

########
# Notes:
# Sample sturctures are associated with the chromosomes.
# Conditional mean structures are usually estimated for each chromosome separately,
# but can later be combined together using sumCondMeans
########

######
# We want to find the best GC window for this data
# We stratify locations based on the %GC of the window starting margin bp after the location.
# We compare different window lengths, in this example between 8 and 400.
#
# (In terms of Benjamini and Speed (2011), we compare W_a,l for a=margin and different l's
# We expect the best value to approximately be the median fragment length
######

# We always recommend throwing away the first few bp's, so as not
# to confound with fragmentation effects. 
margin = 5

begdata = makeGCLens(chr1$isgc,dat1$forw,sampline = sampCh1$singleLocSamp,minlen = 8,maxlens=650, margin=5,max_frag_for_loc=10)

#####
# begdata consists of:
# begdata$locs: a matrix of locations
# begdata$frag: a matrix of fragment counts
# row l in each corresponds W_{0,l} (stratification based on a GC window of size l)
#####

#####
# Computes the TV score for each row (deviation from uniformity)
####
tvs = scoreGCLens(begdata, maxlen=650, minlen =  8,scale= T)
plotGCLens(tvs,lw=2,lt=1)
best_size = which.max(tvs)-1

####
# We now use this size to generate condMean table...
# We want to make a gc line of window-2*margin
####

gcsize = best_size # I got 548 for my sample of the full chromosome
gcline = prepareGCChrom(chr1,gcsize,filename = "gcchr1_sm")


##
# gcline[i] = isgc[i]+isgc[i+1]+....isgc[i+gcsize-1]
# We want to stratify location i based on gcline[i+margin]
# We use our sample for this, comparing the gcline[sampCh1$singleLocSamp+margin]
# to the fragment counts at these locations (already computed in sampCh1$forsamped)
##

cMeans= getCondMean(gcline[sampCh1$singleLocSamp+margin],sampCh1$forsamped,cutoff = 4,jump = 6,norm = FALSE)
# Taking jumps at jumps of 3 give smoother results.

# Similarly, we could estimate the condMean on the reverse strand
# Here, we need to do the adjustment:
# First find the 3' end (x + readlen), then go back -(margin+gcsize)
cMeans_rev= getCondMean(gcline[sampCh1$singleLocSamp+dat1$readlen-1-margin-gcsize]
  ,sampCh1$revsamped,cutoff = 4,jump = 6,norm = FALSE)
  

# And plot
par(mfcol = c(1,2))
plotCondMean(cMean = cMeans,ci = TRUE,normRange = gcsize,index=0,meanLine=TRUE,lt = 1,col=4)
plotCondMean(cMean = cMeans_rev,ci = TRUE,normRange = gcsize,index=1,meanLine=TRUE,lt = 1,col=4)

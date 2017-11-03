module load Python/2.7.9-anaconda 
file="SRR961514"
python coverage.py -1 ${file}_1.fastq -2 ${file}_2.fastq -x indexes/K03455.fna -q 30

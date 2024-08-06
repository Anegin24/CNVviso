#####multiple-mapping-with-bwa#####
#test-on-PE-data-PGS#
#read<70bp#

ref=~/Bioinformatics/Ref/hg38chr.fa

forwardread=(~/Bioinformatics/231222-L65CR/*_R1*)
for forwardread in ${forwardread[@]}
do
reverseread=$(echo $forwardread | sed 's\R1\R2\g')

read1sainame=$(basename $forwardread | cut -f1,2,3,4 -d_)
read2sainame=$(basename $reverseread | cut -f1,2,3,4 -d_)
outputfilename=$(basename $forwardread | cut -f1,2 -d_)

bwa aln $ref $forwardread > ~/Bioinformatics/231222-L65CR/sai/${read1sainame}.sai 
bwa aln $ref $reverseread > ~/Bioinformatics/231222-L65CR/sai/${read2sainame}.sai

read1=~/Bioinformatics/231222-L65CR/sai/${read1sainame}.sai
read2=~/Bioinformatics/231222-L65CR/sai/${read2sainame}.sai

bwa sampe $ref $read1 $read2 $forwardread $reverseread > ~/Bioinformatics/231222-L65CR/sam/${outputfilename}.sam

samfile=~/Bioinformatics/231222-L65CR/sam/${outputfilename}.sam

samfile=~/Bioinformatics/231222-L65CR/sam/${outputfilename}.sam

cd ~/Bioinformatics/231222-L65CR/processbam

samtools view -b -o ${outputfilename}.bam $samfile

samtools view -b -F 0xc ${outputfilename}.bam -o ${outputfilename}.filtered.bam

samtools sort -@ 20 -n ${outputfilename}.filtered.bam -o ${outputfilename}.sorted.n.bam

samtools fixmate -m ${outputfilename}.sorted.n.bam ${outputfilename}.fixmate.bam

samtools sort -@ 20 ${outputfilename}.fixmate.bam -o ${outputfilename}.sorted.p.bam

samtools markdup -r -@ 20 ${outputfilename}.sorted.p.bam ~/Bioinformatics/231222-L65CR/finalbam/${outputfilename}.dedup.bam


done


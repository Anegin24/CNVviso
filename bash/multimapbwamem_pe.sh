#####multiple-mapping-with-bwa#####
#test-on-SE-data-PGS#
#read<70bp#

ref=~/Bioinformatics/Ref/hg38chr.fa

forwardread=(~/Bioinformatics/Test31/fastq/*_R1*)
for forwardread in ${forwardread[@]}
do
reverseread=$(echo $forwardread | sed 's\R1\R2\g')

outputfilename=$(basename $forwardread | cut -f1 -d_)

bwa mem $ref $forwardread $reverseread > ~/Bioinformatics/Test31/sampe/${outputfilename}.sam

cd ~/Bioinformatics/Test31/processbampe

samtools view -b -o ${outputfilename}.bam ~/Bioinformatics/Test31/sampe/${outputfilename}.sam

samtools view -b -F 0xc ${outputfilename}.bam -o ${outputfilename}.filtered.bam

samtools sort -@ 20 -n ${outputfilename}.filtered.bam -o ${outputfilename}.sorted.n.bam

samtools fixmate -m ${outputfilename}.sorted.n.bam ${outputfilename}.fixmate.bam

samtools sort -@ 20 ${outputfilename}.fixmate.bam -o ${outputfilename}.sorted.p.bam

samtools markdup -r -@ 20 ${outputfilename}.sorted.p.bam ~/Bioinformatics/Test31/finalbampe/${outputfilename}.dedup.bam

done

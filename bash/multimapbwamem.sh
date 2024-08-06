#####multiple-mapping-with-bwa#####
#test-on-SE-data-PGS#
#read<70bp#
inputdirectory=~/Bioinformatics/Demo/
ref=~/Bioinformatics/Ref/hg38chr.fa
refcnvkit=~/Bioinformatics/Ref/refcnvkit/hg38.cnn
mkdir $inputdirectory/sam
mkdir $inputdirectory/processbam
mkdir $inputdirectory/finalbam
mkdir $inputdirectory/cnvkitresult

forwardread=($inputdirectory/fastq/*_R1*)
if [ ${#forwardread[@]} -eq 0 ]; then
    echo "No fastq files found. Removing directories."
    rmdir $inputdirectory/sam $inputdirectory/processbam $inputdirectory/finalbam $inputdirectory/cnvkitresult
    exit 1
fi


for forwardread in ${forwardread[@]}
do
#reverseread=$(echo $forwardread | sed 's\R1\R2\g')

outputfilename=$(basename $forwardread | cut -f1 -d_)

cd $inputdirectory

bwa mem $ref $forwardread > $inputdirectory/sam/${outputfilename}.sam

cd $inputdirectory/processbam

samtools view -b -o ${outputfilename}.bam $inputdirectory/sam/${outputfilename}.sam

samtools view -b -F 0xc ${outputfilename}.bam -o ${outputfilename}.filtered.bam

samtools sort -@ 20 -n ${outputfilename}.filtered.bam -o ${outputfilename}.sorted.n.bam

samtools fixmate -m ${outputfilename}.sorted.n.bam ${outputfilename}.fixmate.bam

samtools sort -@ 20 ${outputfilename}.fixmate.bam -o ${outputfilename}.sorted.p.bam

samtools markdup -r -@ 20 ${outputfilename}.sorted.p.bam $inputdirectory/finalbam/${outputfilename}.dedup.bam

rm $inputdirectory/sam/${outputfilename}.sam
rm $inputdirectory/processbam/${outputfilename}*

##########processcnvkit###################

cnvkit.py batch -m wgs $inputdirectory/finalbam/${outputfilename}.dedup.bam -r $refcnvkit -p 8 --scatter --diagram -d $inputdirectory/cnvkitresult/${outputfilename}/

cd $inputdirectory/cnvkitresult/${outputfilename}

zip ${outputfilename}.dedup.zip ${outputfilename}.dedup.cnr ${outpurfilename}.dedup.cns ${outputfilename}.dedup.call.cns ${outputfilename}

done


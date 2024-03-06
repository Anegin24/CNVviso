#!/bin/bash
script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <input_directory>"
    exit 1
fi

# Assign the input directory provided as an argument
inputdirectory="$1"

# Define other paths and directories
ref="$script_dir/Ref/hg38chr.fa"
refcnvkit="$script_dir/Ref/hg38.cnn"

# Create necessary directories
mkdir -p "$inputdirectory/sam"
mkdir -p "$inputdirectory/processbam"
mkdir -p "$inputdirectory/finalbam"
mkdir -p "$inputdirectory/cnvkitresult"

# Retrieve fastq files from the input directory
forwardread=($inputdirectory/fastq/*_R1*)

if [ ${#forwardread[@]} -eq 0 ]; then
    echo "No fastq files found in $inputdirectory. Exiting."
    exit 1
fi

# Iterate through fastq files and perform processing
for forwardread in "${forwardread[@]}"; do
    # Extract the output filename from the fastq file name
    outputfilename=$(basename "$forwardread" | cut -f1 -d_)

    # Run BWA-MEM
    bwa mem "$ref" "$forwardread" > "$inputdirectory/sam/${outputfilename}.sam"

    # Process BAM files
    samtools view -b -o "${outputfilename}.bam" "$inputdirectory/sam/${outputfilename}.sam"
    samtools view -b -F 0xc "${outputfilename}.bam" -o "${outputfilename}.filtered.bam"
    samtools sort -@ 20 -n "${outputfilename}.filtered.bam" -o "${outputfilename}.sorted.n.bam"
    samtools fixmate -m "${outputfilename}.sorted.n.bam" "${outputfilename}.fixmate.bam"
    samtools sort -@ 20 "${outputfilename}.fixmate.bam" -o "${outputfilename}.sorted.p.bam"
    samtools markdup -r -@ 20 "${outputfilename}.sorted.p.bam" "$inputdirectory/finalbam/${outputfilename}.dedup.bam"

    # Clean up intermediate files
    rm "$inputdirectory/sam/${outputfilename}.sam"
    rm "${outputfilename}.bam"
    rm "${outputfilename}.filtered.bam"
    rm "${outputfilename}.sorted.n.bam"
    rm "${outputfilename}.fixmate.bam"

    # Process with CNVkit
    cnvkit.py batch -m wgs "$inputdirectory/finalbam/${outputfilename}.dedup.bam" -r "$refcnvkit" -p 8 --scatter --diagram -d "$inputdirectory/cnvkitresult/${outputfilename}/"
    cd "$inputdirectory/cnvkitresult/${outputfilename}"
    pdftoppm -png -r 200 *.pdf "${outputfilename}_diagram"
    zip "${outputfilename}.dedup.zip" "${outputfilename}.dedup.cnr" "${outputfilename}.dedup.cns" "${outputfilename}.dedup.call.cns" "${outputfilename}_diagram.png" "${outputfilename}"
done

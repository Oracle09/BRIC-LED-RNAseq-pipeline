#!/bin/bash
find . -type f -exec mv {} . \;     #Move all the replicate files from sub directory to main folder

# Run FastQC for quality control
fastqc *.fastq.gz > 'Quality check for all fastq.gz files in folder'

# unzip .gz files which is raw. This is not necessary to save space
for file in *.fastq.gz;
do 
gunzip -k "$file" #-k retains the original file
done

#Trim using Trimgalore single since not paired data


#install python cutadapt before trimming and specify the path
apt install cut adapt


for file in *fastq.gz
do
trim_galore -j 4 -o TrimGalore "$file"   #4 cores should be used   
done

#STAR Alignment 
#download Athaliana_447_Araport11.gene_exons.gff3 as the annotation and Athaliana_447_TAIR10.fa as the reference genome

STAR --runThreadN 4 --runMode genomeGenerate --genomeDir /mnt/c/Users/Gbolaga/Desktop/'BRIC LED RNA RAW DATA'/STAR_index \
--genomeFastaFiles Athaliana_447_TAIR10.fa --sjdbGTFfile Athaliana_447_Araport11.gene_exons.gff3 \ --sjdbOverhang 1
--sjdbOverhang 99

STAR --genomeDir /mnt/c/Users/Gbolaga/Desktop/'BRIC LED RNA RAW DATA'/STAR_index/ --readFilesIn FT-C2_S9_L001_R1_001_trimmed.fq.gz --runThreadN 4 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix /home/Users/Gbolaga/Desktop/'BRIC LED RNA RAW DATA'/STAR_output/
        STAR --genomeDir "/mnt/c/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/STAR_index/" --readFilesIn FT-C2_S9_L001_R1_001_trimmed.fq.gz --runThreadN 4 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix "/home/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/STAR_output/"


#Sort the BAM file by chromosomal position
samtools sort -@ 8 -o /home/Users/Gbolaga/Desktop/BRIC\ LED\ RNA\ RAW\ DATA/STAR_output/Aligned.sortedByCoord.sorted.bam /home/Users/Gbolaga/Desktop/BRIC\ LED\ RNA\ RAW\ DATA/STAR_output/Aligned.sortedByCoord.out.bam
#This retains the .bam file while the newly sorted one contains sorted.out.bam

#Count read the sorted BAM file
featureCounts -a /mnt/c/Users/Gbolaga/Desktop/'BRIC LED RNA RAW DATA'/Athaliana_447_Araport11.gene_exons.gff3 -o counts.txt -g ID /home/Users/Gbolaga/Desktop/BRIC\ LED\ RNA\ RAW\ DATA/STAR_output/Aligned.sortedByCoord.sorted.bam
#specify -g ID because ggf files have ID instead of gene ID


#To exceute this for multiple files in the folder
#!/bin/bash

#Define directories
INPUT_DIR="/mnt/c/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/GE5429-RFerl_S4-H2MHVDSX3-Lane1"
INPUT_DIR2="/mnt/c/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA"
STAR_INDEX_DIR="$INPUT_DIR2/STAR_index"
OUTPUT_DIR="/home/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/STAR_output"
ANNOTATIONS="/mnt/c/Users/Gbolaga/Desktop/BRIC LED RNA RAW DATA/Athaliana_447_Araport11.gene_exons.gff3"

#Loop over .fq.gz files in the INPUT_DIR directory
for fq_file in "$INPUT_DIR"/*_trimmed.fq.gz; do
    #Get the base name for the file
    base_name=$(basename "$fq_file" _trimmed.fq.gz)

    #STAR alignment
    STAR --genomeDir "$STAR_INDEX_DIR" --readFilesIn "$fq_file" --runThreadN 4 --outSAMtype BAM SortedByCoordinate --readFilesCommand zcat --outFileNamePrefix "$OUTPUT_DIR/${base_name}_"

    #Sort BAM file
    samtools sort -@ 8 -o "$OUTPUT_DIR/${base_name}_Aligned.sortedByCoord.out.bam" "$OUTPUT_DIR/${base_name}_Aligned.sortedByCoord.out.bam"

    #Count reads with featureCounts
    featureCounts -a "$ANNOTATIONS" -o "$OUTPUT_DIR/${base_name}_counts.txt" -g ID "$OUTPUT_DIR/${base_name}_Aligned.sortedByCoord.out.bam"
done
Do the average first in R
#To reduce the file names to just 5 characters e.g GC-C2
for file in *; do
  if [[ -f "$file" ]]; then
    extension="${file##*.}"
    shortened_name="${file%%_*}"
    new_name="${shortened_name}.${extension}"
    cp "$file" "$new_name"
  fi
done
#!/bin/bash

# human respiratory and olfactory epithelium
while read sample fastq; do
    echo "Mapping sample ${sample}..."
    STAR --genomeLoad LoadAndKeep \
         --runThreadN 4 \
         --genomeDir /home/irlab/E/star-indexes/Homo_sapiens.GRCh38.99 \
         --readFilesIn ${fastq} \
         --outFilterMultimapNmax 1 \
         --outFileNamePrefix ${sample}.

done <sample_list.txt

# remove the genome from the shared memory
echo "Removing genome from the shared memory"
STAR --genomeLoad Remove --genomeDir /home/irlab/E/star-indexes/Homo_sapiens.GRCh38.99


#!/bin/bash

for FILEIN in *.Aligned.out.sam; do
    featureCounts -a '/home/irlab/E/Annotations/Homo_sapiens.GRCh38.99.no_AC097625.2.gtf' \
                  -o ${FILEIN%Aligned.out.sam}featureCounts_counts.2.txt \
                  -F GTF \
                  -t exon \
                  -g gene_id \
                  --minOverlap 1 \
                  --fracOverlap 0 \
                  --fracOverlapFeature 0 \
                  -Q 0 \
                  -s 2 \
                  -T 1 \
                  -R BAM \
                  -T 4 \
                  $FILEIN;
done

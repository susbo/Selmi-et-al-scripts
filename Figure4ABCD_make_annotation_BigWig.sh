#!/bin/bash

software="/home/USSR/smb208/Projects/Software"
chrom_sizes="/servers/bio-shares/bioinf-facility/smb208/Reference/hg38/Annotation/hg38.chrom.sizes"

wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.basic.annotation.gtf.gz

zcat gencode.v28.basic.annotation.gtf.gz | awk 'OFS="\t" {if ($3=="exon") {print $1,$4-1,$5,$10,$26,$7}}' | tr -d '";' | sort -k 1,1 -k2,2n | uniq > exons.hg38.bed
zcat gencode.v28.basic.annotation.gtf.gz | awk 'OFS="\t" {if ($3=="CDS") {print $1,$4-1,$5,$10,$26,$7}}' | tr -d '";' | sort -k 1,1 -k2,2n | uniq > cdss.hg38.bed

for f in exons cdss
do
	for strand in plus minus
	do
		s="+"
		if [ "$strand" == "minus" ]; then s="-"; fi
		cat $f.hg38.bed | awk -v s=$s '$6==s' > $f.hg38.$strand.bed
		bedtools merge -i $f.hg38.$strand.bed | awk 'BEGIN{OFS="\t"} {print $0,"1"}' > $f.hg38.$strand.merged.bed
   	$software/bedGraphToBigWig $f.hg38.$strand.merged.bed $chrom_sizes $f.hg38.$strand.merged.bw
	done
done

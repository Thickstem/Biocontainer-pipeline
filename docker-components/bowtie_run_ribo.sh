SRRID=SRR1562539
"""
# remove contaminants

bowtie2 --very-sensitive-local -x ./refs/index/contaminant \
		-U ./output/${SRRID}/${SRRID}_trimmed.fastq \
		--un ./output/${SRRID}/${SRRID}_trimmed_filtered.fastq \
		-S ./output/${SRRID}/${SRRID}_contaminant.sam \
		-p 32 
rm ./output/${SRRID}/${SRRID}_contaminant.sam
"""

## main mapping 
bowtie2 --very-sensitive-local -k 100 --norc \
		--rdg 99999999,99999999 --rfg 99999999,99999999 \
		-x ./refs/index/mRNA.lncRNA \
		-U ./output/${SRRID}/${SRRID}_trimmed_filtered.fastq \
		-S ./output/${SRRID}/${SRRID}_main.sam \
		-p 32

samtools view -F4 -S -h ./output/${SRRID}/${SRRID}_main.sam > ./output/${SRRID}/${SRRID}_mapped.sam
samtools view -c ./output/${SRRID}/${SRRID}_mapped.sam

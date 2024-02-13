SRRID=SRR2075925
#ADAPTER1=CTGTAGGCACCAT # for SRR1562539
ADAPTER1=CTGTAGGCACCATCAAT # for SRR2075925
ADAPTER2=AAAAAAAAAAAAAAAAAAAAAAAA

cutadapt -m 25 -q 10 -a ${ADAPTER1} -a ${ADAPTER2} -j 32 -o ./output/${SRRID}/${SRRID}_trimmed.fastq \
		 ./output/${SRRID}/${SRRID}.fastq.gz


## -m : minimum length
## -j : num cores to use

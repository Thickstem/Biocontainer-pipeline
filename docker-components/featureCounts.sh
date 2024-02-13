SRRID=SRR8550317
featureCounts -a ./refs/gencode.v45.annotation.gtf -o ./output/${SRRID}/counts.txt \
				-t transcript -g transcript_id \
				./output/${SRRID}/${SRRID}_ribo.sam

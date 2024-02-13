SRRID=SRR1562539

python read_count.py --fasta ./refs/mRNA.lncRNA.fa \
		--sam ./output/${SRRID}/${SRRID}_ribo_mapped_mod.sam \
		--out_dir ./output/${SRRID}/

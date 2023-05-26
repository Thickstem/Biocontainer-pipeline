IMAGE_NAME=ncbi/sra-tools:3.0.1

run:
	docker run --rm -it --name fastq-test \
				-v ${PWD}:/tmp/wk \
				${IMAGE_NAME} 

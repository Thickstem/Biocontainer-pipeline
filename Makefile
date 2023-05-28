IMAGE_NAME=ncbi/sra-tools:3.0.1

run:
	docker run --rm -it --name fastq-test \
				-v ${PWD}:/tmp/wk \
				${IMAGE_NAME} 

all:
	docker-compose -f docker-components/docker-compose.yml up

sra:
	docker-compose -f docker-components/docker-compose-sratoolkit.yml up

fastq:
	docker-compose -f docker-components/docker-compose-fastq-dump.yml up

fastp:
	docker-compose -f docker-components/docker-compose-fastp.yml up

star:
	docker-compose -f docker-components/docker-compose-star.yml up

star-idx:
	docker-compose -f docker-components/docker-compose-star-index.yml up

rsem:
	docker-compose -f docker-components/docker-compose-rsem.yml up

rsem-ref:
	docker-compose -f docker-components/docker-compose-rsem-ref.yml up




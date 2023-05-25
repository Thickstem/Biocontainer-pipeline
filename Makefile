IMAGE_NAME=biocontainers/rsem

run:
	docker run --rm -it --name rsem-test \
				-v ${PWD}:/tmp/wk \
				${IMAGE_NAME} bash

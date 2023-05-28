# BioContainer pipeline
- Build a pipeline create transcript level count matrix from SRR file only with Biocontainers
- This pipeline can be used any environment containing docker/docker-compse. (Win/Mac/Linux)
## Dir structure
```
Biocontainers
|--Dockerfiles
|--input
|--output
|--STARidx
|--docker-components
|		|- docker-compose-XXX.yml
|		|-.env
|--wait-for-it.sh

```
## Overall flow
- All processes are cotroled in `docker-compose.yml`
- `wait-for-it.sh` dominates the run order of containers written in  docker-compose 
```
SRR file
↓ 1.(parallel-fastq-dump)
FastQ 
↓ 2.(fastp)
FastQ (quarity checked.trimed)
↓ 3.(STAR)
BAM/SAM
↓ 4.(RSEM)
count data
```

## Usage
1. create docker image for RSEM
```sh
docker build -t biocontainer/rsem -f Dockerfiles/rsem
```
2. Create STAR index (refer to [this section](https://github.com/Thickstem/Biocontainer-pipeline/blob/master/Details.md#31-prepare-index) )
3. Create RSEM index (refer to [this section](https://github.com/Thickstem/Biocontainer-pipeline/blob/master/Details.md#41-prepare-reference))
4. Modify variables in `.env` 

### All run
```sh
make all
```
### Each tool run
- This repository has `docker-compose-XXX.yml` file for each bioinformatics tools
	- XXX: sratoolkit,fastq-dump,fastp,STAR,RSEM
- Each tools can be run independently with below command
```sh
make ${TOOL_NAME}
```

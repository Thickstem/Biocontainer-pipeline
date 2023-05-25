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
|--.env
|--docker-compose.yml
|--wait-for-it.sh

```
## Usage
1. create docker image for RSEM
```sh
docker build -t biocontainer/rsem -f Dockerfiles/rsem
```
2. Create STAR index
3. Create RSEM index


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


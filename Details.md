## 1.fastq-dump
- Create FASTQ file from SRR file.
- SRR file have to be dowloaded and put to `./input` before run.
	- SRR can't download directry with fastq-dump because of permission deny to NCBI.
- This command is for single-end. Modify options when running for pair-end.
### Command
```sh
parallel-fastq-dump --sra-id ./input/${ID} --threads ${THREADS} --gzip &&
      mkdir -p output/${ID} &&
      mv *.fastq.gz output/${ID}/"
```
- `--sra-id`:SRA file path
- `--thread`: Num of thread to use (max:32 @leap-bis)
- `--gzip`: it generate .fastq.gz file

## 2.FastP
- It preprocess for FASTQ file.
    - Triming of Adapter sequence
    - Filtering low quality reads
### Command
```sh
fastp -w 16 -h output/${ID}/${ID}.html -j output/${ID}/${ID}.json -i output/${ID}/${ID}.fastq.gz  -o output/${ID}/${ID}_trim.fastq.gz 
```
- `-w`:num workers (max 16. this maximum is depends on fastp software itself)
- `-h`:Save path of html format result
- `-j`:Save path of json format result
- `-i`:Input file path. Generated from 1.fastq-dump
- `-o`:Outpuf file path. It should contain XX_trim

## 3.STAR
### 3.1 Prepare Index
- STAR need to prepare Index before mapping run.
- This process needs two file. It is recommended to download that two file from [GENCODE](https://www.gencodegenes.org/human/)
	- GTF/GFF3: [Comprehensive gene annotation (PRI)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/gencode.v43.primary_assembly.annotation.gtf.gz)
	- FASTA: [Genome sequence, primary assembly (GRCh38) (PRI)](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_43/GRCh38.primary_assembly.genome.fa.gz)
* When considering only mRNA mapping, Primary region is better because of file size.
#### Command
```sh
STAR --runThreadN ${THREADS} \
 	 --runMode genomeGenerate \
	 --genomeDir STARidx \
	 --genomeFastaFiles STARidx/${FASTA file name} \
	 --sjdbGTFfile STARidx/${GTF/GFF3 file name} \
	 ${STARidxopt}
```

### 3.2 Mapping
#### Command
```sh
STAR --runThreadN ${THREADS} \
	--genomeDir STARidx \
	--readFilesIn output/${ID}/${ID}_trim.fastq.gz \
	--readFilesCommand gunzip -c \
	--quantMode TranscriptomeSAM \
	--genomeLoad NoSharedMemory \
	--outSAMtype BAM SortedByCoordinate \
	--outFileNamePrefix output/${ID}/${ID}_trim. \
	--outFilterMultimapNmax 1
```
- If you want to use RSEM, `--quantMode TranscriptomeSAM` is necessary.

## 4.RSEM
- It count read mapping result transcript level.
### 4.1 Prepare-reference
- RSEM needs reference files (different from STAR or other index)
- Run count 
#### Command
```sh
rsem-prepare-reference --gtf ${GTF file path} \
					   -p ${THREADS} \
					   FASTA_file_path \
					   ref_name_prefix
```
- `--gtf`: Path to GTF file. `--gff3` is also able.
- `-p`: process num. = Threads num 
- `FASTA_file_path`: .fa file path. It should be the same file used in STAR index preparation.
- `ref_name`: RSEM reference output some ref files with this predix. 

### 4.2 Count
```sh
rsem-calculate-expression \
          --alignments \
          --append-names \
          -p ${THREADS} \
          Path/to/XXX.toTranscriptome.out.bam \
          ref_name_prefix \
          sample_name
```

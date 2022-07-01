# RNA Alternative Splicing Analyses

This page is the collection of how I learn and conduct analyses for bulk RNA sequence. 

The main analyses here, is to use the raw sequence for mapping, and then do the **Differential Gene Expression analysis**, **Alternative Splicing Analysis**, **Motif Enrichment Analysis** and **Chromosome Analysis**.


## Sequence Preparation: Kallisto 

Before doing the **Differential Gene Expression Analysis**, we need to use Kallisto to do the bulk RNA sequence mapping to get the raw counts of events.


### First, download reference genome files from Ensembl webpage both all.fa.gz file and gtf.gz file in to the Kallisto directory and then extract them.

Ensembl webpage: https://useast.ensembl.org/info/data/ftp/index.html


Command:
`wget -P /path/where/you/want/the/file/to/be/downloaded/ http://~(file link)`

For example ( Human gene ):

*Get .fa.gz*.        
`wget -P /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/ http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`

*Extract to .fa*.     
`gunzip /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.cdna.all.fa.gz`

*Get gtf.gz*.  
`wget -P /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/ http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz`

*Extract to .gtf*.  
`gunzip /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.104.gtf.gz`

### Second, generate transcriptome.idx in Kallisto directory

Using batch file: **KallistoTranscriptome.sh**

Command:
`sbatch /path/to/the/batch/file/KallistoTranscriptome.sh`

```markdown
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00 # 24 hr
#SBATCH --mem=40g
#SBATCH -c 4 #4 cpus

module load kallisto

kallisto index -i /where/you/want/to/store/Kallisto/[Hsapiens or Mmusculus]/transcriptome.idx /reference/file/downloaded/in the/first/step/[filename].cdna.all.fa

# Example
kallisto index -i /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/transcriptome.idx /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.cdna.all.fa

```
### Third, download sample reads file (.fastq.gz)

Get the sequence fastq files

### Fourth, take RNA-sequencing samples (.fastq.gz files) and align them to the Kallisto reference

Using batch file: **Kallisto.sh**

Command:
`sbatch /path/to/the/batch/file/Kallisto.sh`

```markdown
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00 # 24 hr
#SBATCH --mem=40g
#SBATCH -c 4 #4 cpus

module load kallisto

kallisto quant -i /where the/reference/file/Kallisto/[Hsapiens or Mmusculus]/transcriptome.idx --pseudobam -b 50 -t 4 -o /where/you/want/to store/ /sample read path/[reads1]_R1_001.fastq.gz /sample read path/[reads2]_R2_001.fastq.gz

# Example
kallisto quant -i /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/transcriptome.idx --pseudobam -b 50 -t 4 -o /home/lk627/project/Test/Lab/Project/Sample#1/ /home/lk627/project/Test/Lab/SampleReads/H2291_1_261_316_S25_L001_R1_001.fastq.gz /home/lk627/project/Test/Lab/SampleReads/H2291_1_261_316_S25_L001_R2_001.fastq.gz

```




### Differential Gene Expression Analysis

Markdown is a lightweight and easy-to-use syntax for styling your writing. It includes conventions for

```markdown
Syntax highlighted code block

# Header 1
## Header 2
### Header 3

- Bulleted
- List

1. Numbered
2. List

**Bold** and _Italic_ and `Code` text

[Link](url) and ![Image](src)
```



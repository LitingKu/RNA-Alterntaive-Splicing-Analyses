# RNA Alternative Splicing Analyses

This page is the collection of how I learn and conduct analyses for bulk RNA sequence. 

The main analyses here, is to use the raw sequence for mapping, and then do the **Differential Gene Expression analysis**, **Alternative Splicing Analysis**, **Motif Enrichment Analysis** and **Chromosome Analysis**.


## Sequence Preparation: Kallisto 

### First, download reference genome files from Ensembl webpage both all.fa.gz file and gtf.gz file in to the Kallisto directory and then extract them.

Ensembl webpage: https://useast.ensembl.org/info/data/ftp/index.html


Command:
`wget -P /path/where/you/want/the/file/to/be/downloaded/ http://~(file link)`

For example ( Human gene ):

**Get .fa.gz**.       
`wget -P /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/ http://ftp.ensembl.org/pub/release-104/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz`

**Extract to .fa**.     
`gunzip /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.cdna.all.fa.gz`

**Get gtf.gz**.  
`wget -P /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/ http://ftp.ensembl.org/pub/release-104/gtf/homo_sapiens/Homo_sapiens.GRCh38.104.gtf.gz`

**Extract to .gtf**.  
`gunzip /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.104.gtf.gz`

### Second, generate transcriptome.idx in Kallisto directory

Using batch file: **KallistoTranscriptome.sh**

```markdown
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00 # 24 hr
#SBATCH --mem=40g
#SBATCH -c 4 #4 cpus

module load kallisto
kallisto index -i /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/transcriptome.idx /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/Homo_sapiens.GRCh38.cdna.all.fa

```



1. Before doing the **Differential Gene Expression Analysis**, I use Kallisto to do the bulk RNA sequence mapping to get the raw counts of events.

`module load kallisto`

`kallisto quant -i /home/lk627/project/Test/Lab/Reference/Kallisto/Hsapiens/transcriptome.idx --pseudobam -b 50 -t 4 -o /gpfs/ycga/scratch60/escobar-hoyos/lk627/SF3B1/Mia/MiaEV/ /gpfs/ycga/scratch60/escobar-hoyos/lk627/SF3B1/Mia/fastq_file/MiaEV1_IGO_10215_B_11_S26_R1_001.fastq.gz /gpfs/ycga/scratch60/escobar-hoyos/lk627/SF3B1/Mia/fastq_file/MiaEV1_IGO_10215_B_11_S26_R2_001.fastq.gz`

```markdown

`module load kallisto`

## EV


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

For more details see [Basic writing and formatting syntax](https://docs.github.com/en/github/writing-on-github/getting-started-with-writing-and-formatting-on-github/basic-writing-and-formatting-syntax).

### Jekyll Themes

Your Pages site will use the layout and styles from the Jekyll theme you have selected in your [repository settings](https://github.com/LitingKu/RNA-Alterntaive-Splicing-Analyses/settings/pages). The name of this theme is saved in the Jekyll `_config.yml` configuration file.

### Support or Contact

Having trouble with Pages? Check out our [documentation](https://docs.github.com/categories/github-pages-basics/) or [contact support](https://support.github.com/contact) and we’ll help you sort it out.

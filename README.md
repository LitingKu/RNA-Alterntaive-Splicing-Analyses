[Home](https://litingku.github.io) |   | [💻 Lab Work](https://litingku.github.io/about.html)

# RNA Alternative Splicing Analyses

This page is the collection of how I learn and conduct analyses for bulk RNA sequence. 

Here shows the preparation for RNA sequence analyses

The main analyses, is to use the raw sequence for mapping, and then do the **Differential Gene Expression analysis**, **Alternative Splicing Analysis**, **Motif Enrichment Analysis** and **Chromosome Analysis**.


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

## Alternative Splicing Detection: SUPPA2

SUPPA2 help us to find the alternative splicing events with the output file of Kallisto(the abundance.tsv file).

Go to the https://github.com/comprna/SUPPA/releases/tag/v2.3 download Source.code(tar.gz) on local computer or uploaded to Open OnDemand /home/Your netid/(HPC cluster). Next, run command: tar -xvzf /path/to/SUPPA-2.3.tar.gz to extract files. So, your should have a folder called SUPPA-2.3 and now suppy.py is in the /SUPPA-2.3/ directory!

In order to successfully run suppy.py, some packages for python need to be installed first:

`pip install panadas`
`pip install skilearn`
`pip install statsmodels`

### Create SUPPA2 table for downstream analyses

Upload *SUPPA_Complier.py* and *Kallisto2Suppa.R* to the directory.

Remember to change the path name and the name for the samples in the *SUPPA_Complier.py* file.
```markdown
## Code need to be modified in SUPPA_Complier.py

basepath = "/**change to the path where the samples located**/"
data_sample = pd.DataFrame()
for sample in [name for name in os.listdir(basepath) if os.path.isdir(os.path.join(basepath, name)) and name.startswith('**change the character**')]:
```

Using batch file: **SUPPA2.sh**

Command:
`sbatch /path/to/the/batch/file/SUPPA2.sh`

The batch file SUPPA2.sh includes codes for:

1. Convert Kallisto to SUPPA TPM inputs
2. Obtain the Gencode reference genome
3. Calculate PSI values per event
4. Generate SUPPA_merge.csv

Result files after running  **SUPPA2.sh**

1. Exists abundance.txt in /your/sample/folder/
2. Exists lots of .ioe and .gtf files in /Reference/SUPPA2/Hsapiens or Mmusculus/
3. Exists A3 A5 AF AL MX RI SE.psi in /your/sample/folder/
4. Exists SUPPA_merge.csv in /your/directory/contains every sample folder/

## Sequecne Alignment: STAR

With using STAR: https://github.com/alexdobin/STAR for sequence alignment, it help us to generate the **.bam** file that we can do downstream analyses for finding novel alternative splicing events by using rMATs or outrigger.

### First, Download both primary assembly (PRI) version of the genome gtf.gz and fa.gz files for in STAR directory

Human: https://www.gencodegenes.org/human/

Mouse: https://www.gencodegenes.org/mouse/


Command:
`wget -P /path/where/you/want/the/file/to/be/downloaded/ http://~(file link)`

For example ( Human gene ):

*Get gtf.gz*.    
`wget -P /home/lk627/project/Test/Lab/Reference/STAR/Hsapiens/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz`

*Extract to .gtf*.     
`gunzip /home/lk627/project/Test/Lab/Reference/STAR/Hsapiens/gencode.v38.primary_assembly.annotation.gtf.gz`

*Get gtf.gz*.      
`wget -P /home/lk627/project/Test/Lab/Reference/STAR/Hsapiens/ http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/GRCh38.primary_assembly.genome.fa.gz`

*Extract to .gtf*.    
`gunzip /home/lk627/project/Test/Lab/Reference/STAR/Hsapiens/GRCh38.primary_assembly.genome.fa.gz`

### Second, generate the STAR index files

Using batch file: **STAR_index.sh** 

Command:

`sbatch /path/to/the/batch/file/STAR_index.sh`

```markdown
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00 # 24 hr
#SBATCH --mem=40g
#SBATCH -c 4 #4 cpus

## module load is to be used on cluster
module load STAR/2.7.9a-GCCcore-10.2.0

STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /where you want to put the index/ --genomeFastaFiles /the loaction of/[ name].primary_assembly.genome.fa --sjdbGTFfile /the location of/[name].primary_assembly.annotation.gtf --sjdbOverhang 100 --genomeChrBinNbits 18 --genomeSAindexNbases 13 --genomeSAsparseD 3

# Example
STAR --runThreadN 8 --runMode genomeGenerate --genomeDir /gpfs/ycga/scratch60/escobar-hoyos/lk627/Mmusculus --genomeFastaFiles /gpfs/ycga/scratch60/escobar-hoyos/lk627/Mmusculus/GRCm38.primary_assembly.genome.fa --sjdbGTFfile /gpfs/ycga/scratch60/escobar-hoyos/lk627/Mmusculus/gencode.vM10.primary_assembly.annotation.gtf --sjdbOverhang 100 --genomeChrBinNbits 18 --genomeSAindexNbases 13 --genomeSAsparseD 3
```

### Third, generate bam file

Using batch file: **STAR.sh** 

Command:

`sbatch /path/to/the/batch/file/STAR.sh`

```markdown
#!/bin/bash
#SBATCH --mail-type=ALL
#SBATCH -t 24:00:00 # 24 hr
#SBATCH --mem=40g
#SBATCH -c 4 #4 cpus

## module load is to be used on cluster
module load STAR/2.7.9a-GCCcore-10.2.0 


STAR --genomeDir /STAR reference location/ \
--runThreadN 8 \
--readFilesIn <(gunzip -c  /fastq file location/reads1.fastq.gz) <(gunzip -c /fastq file location/reads2.fastq.gz) \
--outFileNamePrefix /the location of the output bam file/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outSAMstrandField intronMotif


# Example
STAR --genomeDir /home/lk627/project/Test/Lab/Reference/STAR/Hsapiens/ \
--runThreadN 8 \
--readFilesIn <(gunzip -c  /home/lk627/project/SRR12900748/SRR12900748_dbGaP-30506_pass_1.fastq.gz) <(gunzip -c / /home/lk627/project/SRR12900748/SRR12900748_dbGaP-30506_pass_2.fastq.gz) \
--outFileNamePrefix /home/lk627/project/SRR12900748/ \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outSAMattributes Standard \
--outSAMstrandField intronMotif
```

## Now all the preparation is done!!

To sum up, with the output of Kallisto we can do the **Differential Gene Expression Analysis** with DESeq2, see the **DESeq2.rmd(DESeq2.html)**

With the output of SUPPA2, we can do **Alternative Splicing Analysis** and **Motif Enrichment Analysis**.

With the output of STAR alignment, we can do sashimi plots with [ggsashimi](https://github.com/guigolab/ggsashimi) and detect novel events by [rMATs](https://github.com/Xinglab/rmats-turbo/blob/v4.1.2/README.md) and [outrigger](https://yeolab.github.io/outrigger/outrigger.index.html#module-outrigger.index).









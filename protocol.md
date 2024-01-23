# Microbial Omics Course (biol217) Protocol - Nils Groeters

# Day 1: Introduction to Linux, Markdown, Git and Conda

### 1 Linux

What we have learned so far:

1. Basic Linux
2. Bioinformatic basics
3. Linux commands

- copy from one folder to another:

```sh
cp source destination
```

# Day 2: Introduction to Metagenomics and Metagenomic Workflow

### Metagenomics

Only a small percentage of the global microbial diversity can be cultured in the lab. In order to study microbes in their natural environment, metagenomic analysis can be used. By sequencing of the genomic DNA, plasmids or transcriptomes it is possible to investigate the microbial diversity.


- Limitations: only proportions, not absolut values of microbes in the environment (see slide 16)

Amplicon Sequencing

- 16S ribisome: 2 gene fused together, very conserved over many domains with conserved and "hyper-variable" regions, conserved regions can be targeted by primers to amplify the vaaible regions
- merenlab.org/momics
- Limitation: only known regions betwenn primers are amplified, mostly 16S RNA (missing plasmids, phages, eucaryotes)
--> shotgun-metagenomic sequencing captures alle sequences

Metagenomic Workflow
(slide 40, 41)

 Sequencing:

- quality control (slide43)
- Phred score for probability of incorrect sequncing result (>20 is good)
- **fastq** file format = read file: contains short reads (up to 150bp)
- sequncing trimming for longer reads with bad quality at the end (phasing or signal decay)

Assembly:

- refernce-guided: known sequences from reference
- de-novo: overlap or de Bruijn

## Assembly of Short Read Sequences on the CAU Cluster (Workflow)

Connect to the CAU cluster by ssh and change to the working directory:

```bash
ssh -X sunam230@caucluster.rz.uni-kiel.de

cd $WORK
```

*Note: For every job a bash script needs to be submitted according to the CAU cluster documentation ([CAU-Cluster](https://www.hiperf.rz.uni-kiel.de/caucluster/)). In order to use the following prigrams the anvio-8 environment needs to be loaded by conda every time a new job is submitted.*

### Program Documentations

| Tool | Version | Repository |
| --- | --- | --- |
| fastqc | 0.12.1 | [FastQC](https://github.com/s-andrews/FastQC ) |
| fastp | 0.23.4 | [fastp](https://github.com/OpenGene/fastp ) |
| megahit | 1.2.9 | [megahit](https://github.com/voutcn/megahit ) |
| samtools | 1.19 | [samtools](https://github.com/samtools/samtools ) |
| QUAST | 5.2.0 | [quast](https://quast.sourceforge.net/quast ) |
| Bowtie2 | 2.4.5 | [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml ) |
| binsanity | 0.5.3 | [binsanity](https://github.com/edgraham/BinSanity) |
| MetaBAT2 | 2.12.1 | [Metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/ ) |
| DASTool | 1.1.5 | [DAS_Tool](https://github.com/cmks/DAS_Tool ) |
| anvi´o | 8 | [anvi’o](https://anvio.org/ ) |
| GUNC | 1.0.5 | [GUNC](https://grp-bork.embl-community.io/gunc/ ) |


**1.** First task is to run multiple raw read files with the **fastqc** program. This program analysis the short read sequences and visualizes the quality scores. The command is included inside the bash script, containing the job, and is handed in on the cau cluster.

Open a bash script file with the emacs editor:

```bash
emacs [filename] &
```

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastqc
#SBATCH --output=fastqc.out
#SBATCH --error=fastqc.err
#SBATCH --partition=all
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam230/[SOURCE DIRECTORY]

# a for loop is used to loop over i files with the ending *fastq.gz
for i in *fastq.gz; do fastqc $i -o [OUTPUT FILE]; done
```
Submit the job:

```bash
sbatch [filename]
```

Tracking the job on the cluster:
 ```bash
 squeue -u [USERNAME]
 ```

**2.** The next task is to run the **fastp** program in order to clean up the raw read files. This is done in pairs (R1/R2), which result from the sequencing method in both directions. **fastp** can not run in a loop over multiple files and has to be executed seperatly for each file pair. The information on the quality of the reads (obtained from **fastqc**) is used for the following setting.

Command line for fastp:
```
fastp -i ? -I ? -R ? -o ? -O ? -t 6 -q 20

> `--html` creates an .html report file in html format\
>`-i` R1 input file name\
>`-I` R2 input file name\
>`-R` report title, here ‘_report’ is added to each file\
>`-o` output_folder/R1.fastq.gz output file\
>`-O` output_folder/R2.fastq.gz output file\
>`-t` trim tail 1, default is 0, here 6 bases are trimmed\
>`-q` 20 reads with a phred score of <=20 are trimmed
```

```bash
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10G
#SBATCH --time=5:00:00
#SBATCH --job-name=fastp
#SBATCH --output=fastp.out
#SBATCH --error=fastp.err
#SBATCH --partition=base
#SBATCH --reservation=biol217

module load gcc12-env/12.1.0
module load miniconda3/4.12.0
conda activate anvio-8

cd /work_beegfs/sunam230/[SOURCE DIRECTORY]
fastp -i BGR_130708_mapped_R1.fastq.gz -I BGR_130708_mapped_R2.fastq.gz -R _report -o ../2_fastp/BGR_130708_mapped_clean_R1.fastq.gz -O ../2_fastp/BGR_130708_mapped_clean_R2.fastq.gz -t 6 -q 20

fastp -i [R1 FILE] -I [R2 FILE] -o [CLEANED R1 FILE] -O [CLEANED R2 FILE] -t -q
...
...
```

**3.** The assembly of the cleaned reads is done by **megahit**. The program takes in all R1 and R2 short read pairs and the output is directed to one folder:

```bash
megahit -1 [R1 FILE] -1 [R1 FILE] -1 [R1 FILE] -2 [R2 FILE] -2 [R2 FILE] -2 [R2 FILE] --min-contig-len 1000 --presets meta-large -m 0.85 -o [OUTPUT FOLDER] -t 12        
```


# Adriana Esparza Pipeline
# This repo will be designated to my LUC COMP BIO python pipeline project- Track 2 Genome Assembly (SPRING 2024)

# Background: For this project, we want to compare HCMV transcriptomes 2 and 6 days post-infection.

# PART 1
# For step 1 I retrieved the four HCMV transcriptomes from two patient donors from SRA and converted them to paired-end fastq files. For each transcriptome I went to NCBI by following each link provided, clciked on the Run ID to get the data, clicked on the Data access tab, then copied link address, and in the terminal used the wget function to retrieve this SRA normalized data via the command line. All transcriptomes I downloaded in a folder within my directory in the class server. To uncompress the data, I used the fastq-dump -I --split-files command which reads each transcriptome file and splits paired-end reads in two files. This last step takes a few minutes to run and produces 2 paired-end files per transcriptome. I then shortened the files to make samples of them and smooth run time when running the code. For this, I used, head -n 40000 data.fastq > nameofmyfile.fastq (and repeated for each of my 8 .fastq files).

# For this project, I chose TRACK 2- Genome Assembly.
# PART 2
# Before assembly, the goal is to map the reads to the HCMV genome. So I first created an index. To do this, I went to NCBI and looked for the genome that matched the entry NC_006273.2. Here is the link for it, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000845245.1/ . When on the website, I clicked on genome, then datasets. This opens a tab with the Datasets command-line query. This can be copied and pasted it to the terminal as follows,
datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report
# This will download the ncbi dataset zip file and validate it. Now it can be unzipped using the following command,
unzip ncbi_dataset.zip
# Now, the index can be created. This will be necessary before mapping. The following command uses Bowtie2-build, then the dataset must be written as well as the name of the index (it can be any but I named mine HCMVINDEX)
bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMVINDEX

#STOPPED HERE.....

# Now we are ready to map and we can use the following command. Note that the file names will have to be changed accordingly.
bowtie2 --quiet -x HCMVINDEX -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMVmap4.sam
# To save only the reads that map,
bowtie2 -x HCMVINDEX -1 SRR5660030_1.fastq -2 SRR5660030
_2.fastq --al-conc mapped_reads.fastq --un-conc unmapped_reads.fastq > PipelineProject.log
# To output the number of reads that each dnor had before and after filtering,
import os



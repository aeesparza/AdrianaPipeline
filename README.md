# AdrianaPipeline
# This repo will be designated to my COMP BIO python pipeline project (SPRING 2024)
# For step 1 I retrieved the four HCMV transcriptomes from two patient donors from SRA and converted them to paired-end fastq files. For each transcriptome I went to NCBI by following each link provided, clciked on the Run ID to get the data, clicked on the Data access tab, then copied link address, and in the terminal used the wget function to retrieve this SRA normalized data via the command line. All transcriptomes I downloaded in a folder within my directory in the class server. To uncompress the data, I used the fastq-dump -I --split-files command which reads each transcriptome file and splits paired-end reads in two files. This last step takes a few minutes to run.

# For this project, I chose TRACK 2- Genome Assembly.
# I first went to NCBI and looked for the genome that matched the entry NC_006273.2. Here is the link for it, https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000845245.1/ . When on the website, I clicked on genome, then datasets. This opens a tab with the Datasets command-line query. Copy and paste it to the terminal as follows,
datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report
# This will download the ncbi dataset zip file and validate it. Now it can be unzipped using the following command,
unzip ncbi_dataset.zip
# Now, we can create an index, which will be necessary before mapping. We can use the following command (this uses Bowtie2-build, then there is the dataset
bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMVINDEX
# Now we are ready to map and we can use the following command. Note that the file names will have to be changed accordingly.
bowtie2 --quiet -x HCMVINDEX -1 SRR5660045_1.fastq -2 SRR5660045_2.fastq -S HCMVmap4.sam

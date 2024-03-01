# Adriana Esparza Pipeline
# This repo will be designated to my LUC COMP BIO python pipeline project- Track 2 Genome Assembly (SPRING 2024)
# Dependencies include, Python and standard libraries like sys, os, and subprocess, Bowtie2, SPAdes, and BLAST+.
# Files needed in the same directory or environment (all attached to this repo and less than 1MB each): Samples of the donors pre-prepped and paired-end.fastq -8 in total, multiple files that compose the index of the NC_006273.2 entry, multiple files that compose the local database of the HCMV subfamily. The python wrapper was also attached as a .py file to download if desired but written below as part of the pipeline.
# I have added a README file to this repo and used it to push my script and all files from the class server to GitHub. To access it and clone it locally, the following command can be used in the terminal: 
git clone https://github.com/aeesparza/AdrianaPipeline.git

# Background: For this project, we want to compare HCMV transcriptomes from 2 donors 2- and 6- days post-infection. 4 transcriptomes in total.

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

# For part 5, I created a local database limited to members of the betaherpesvirinae subfamily. For this, I first searched for it in NCBI, pre-selected Refseq to filter the results, and when found clicked on Send to, and then selected Fasta for file type. This will save the file to your computer (I named my file Hsubfamily.fasta) so I then saved it to my working directory. Once there, I used the following command to create the local database which I named betaherpesvirinae. In the terminal:
makeblastdb -in Hsubfamily.fasta -out betaherpesvirinae -title betaherpesvirinae -dbtype nucl

# Continued PART 2-5
# Now we can move on to the python wrapper!

# ADRIANA ESPARZA LUC COMPBIO 483 
# This python wrapper script takes samples of 2 donors and a total of 4 transcriptomes.
# The command-line used was the following (but may vary based on sample input data): python AE_pywrapper.py HCMVINDEX python PANIC.py HCMVINDEX 30sample_reads_1.fastq 30sample_reads_2.fastq 33sample_reads_1.fastq 33sample_reads_2.fastq 44sample_reads_1.fastq 44sample_reads_2.fastq 45sample_reads_1.fastq 45sample_reads_2.fastq
# The code expects the command-line to provide the script's name, a previously created index (in this case for HCMV NCBI accession NC_006273.2), and 8 paired-end files (based on the 4 transcriptomes split in fastq format. Each is less than 1MB)
# Error comments and exits were added just to identify issues while building this script

# ADRIANA ESPARZA LUC COMPBIO 483 
# This python wrapper script takes samples of 2 donors and a total of 4 transcriptomes.
# The command-line used was the following (but may vary based on sample input data): python AE_pywrapperFinal.py HCMVINDEX 30sample_reads_1.fastq 30sample_reads_2.fastq 33sample_reads_1.fastq 33sample_reads_2.fastq 44sample_reads_1.fastq 44sample_reads_2.fastq 45sample_reads_1.fastq 45sample_reads_2.fastq
# The code expects the command-line to provide the script's name, a previously created index (in this case for HCMV NCBI accession NC_006273.2), and 8 paired-end files (based on the 4 transcriptomes split in fastq format. Each is less than 1MB)
# Error comments and exits were added just to identify issues while building this script

import sys
import os
import subprocess

directory_name = "PipelineProject_Adriana_Esparza"

# To create my directory if it doesn't exist (but should exist already)
if not os.path.exists(directory_name):
    os.mkdir(directory_name)

# To move within my directory
os.chdir(directory_name)

# To solve PART 2, MAPPING TO INDEX WITH BOWTIE2 AND COUNTING BEFORE/AFTER READS
def count_aligned_reads(sam_file, mapq_threshold=20): # To count the reads in sam files above a mapq score of 20 *for a 99% confidence level
    aligned_reads = 0
    with open(sam_file, 'r') as f:
        for line in f:
            if line.startswith('@'):
                continue  # To skip header lines containing @
            parts = line.split('\t')
            mapq = int(parts[4])  # Since MAPQ is the fifth field in SAM format this will detect it
            if mapq >= mapq_threshold:
                aligned_reads += 1
    return aligned_reads

def run_bowtie2_and_count_reads(read_file_1, read_file_2, index_prefix): # To calculate the reads before alignment from a FASTQ file (counts lines and divides by 4)
    total_reads_before = sum(1 for line in open(read_file_1)) // 4
    total_reads_before += sum(1 for line in open(read_file_2)) // 4

    sam_output = f"{os.path.splitext(read_file_1)[0]}_bowtie2.sam" # To create a command to align the reads to the index. It then saves the outputs in .sam files
    unaligned_output_base = os.path.splitext(read_file_1)[0] + "_unaligned" # To separate the unaligned reads in a different file

    # To run Bowtie2
    bowtie2_cmd = [
        "bowtie2",
        "-x", index_prefix, # To identify my index and use it to search for alignment locations
        "-1", read_file_1, # To specify that the sample files provided are paired-end
        "-2", read_file_2, # Same here but for the .2s files
        "--un-conc", f"{unaligned_output_base}%.fastq", # To separate the unaligned reads from the aligned ones
        "-S", sam_output # o output the files in the regularly used sam format
    ]

    try:
        subprocess.run(bowtie2_cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error running Bowtie2: {e}")
        sys.exit(1)  # To make the code exit if Bowtie2 fails

    aligned_reads = count_aligned_reads(sam_output, mapq_threshold=20) # To count the aligned reads from the .sam files

    return total_reads_before, aligned_reads

def process_paired_end_reads(read_files, index_prefix): # This performs the before/after filtering counting
    read_pairs_before = {} # The counting is performed in pairs since my files come from only two donors
    read_pairs_after = {}
    for i in range(0, len(read_files), 2):
        donor_id = (i // 2) // 2 + 1  # This will just adjust per donor 1 or 2
        key = f"Donor {donor_id}"
        before, after = run_bowtie2_and_count_reads(read_files[i], read_files[i + 1], index_prefix)
        if key in read_pairs_before:
            read_pairs_before[key] += before
            read_pairs_after[key] += after
        else:
            read_pairs_before[key] = before
            read_pairs_after[key] = after
    return read_pairs_before, read_pairs_after

def write_to_log(filename, read_pairs_before, read_pairs_after): # Later, the output log file will be defined so this will just make sure to write the answer as instructed in there
    with open(filename, "w") as log_file:
        for key in sorted(read_pairs_before.keys()):
            log_file.write(f"{key} had {read_pairs_before[key]} read pairs before Bowtie2 filtering and {read_pairs_after[key]} read pairs after.\n")

# To convert the produced .sam files to .bam and then to .fastq
def sam_to_fastq(sam_file, fastq_file):
    try:
        with open(sam_file, "r") as sam, open(fastq_file, "w") as fastq:
            for line in sam:
                if line.startswith("@"):
                    continue  # Will skip the header lines
                fields = line.split("\t")
                read_id = fields[0]
                sequence = fields[9]
                quality = fields[10]
                fastq.write(f"@{read_id}\n{sequence}\n+\n{quality}\n")
    except Exception as e:
        print(f"Error converting {sam_file} to FASTQ: {e}")
        return False
    return True

def convert_sam_files_to_fastq(sam_files):
    for sam_file in sam_files:
        fastq_file = os.path.splitext(sam_file)[0] + '.fastq'
        if not sam_to_fastq(sam_file, fastq_file):
            print(f"Conversion failed for {sam_file}.")
            return False
        print(f"Conversion successful. FASTQ file saved as {fastq_file}")
    return True

sam_files = ['30sample_reads_1_bowtie2.sam', '33sample_reads_1_bowtie2.sam', '44sample_reads_1_bowtie2.sam', '45sample_reads_1_bowtie2.sam']
convert_sam_files_to_fastq(sam_files)

# To solve PART 3, ASSEMBLING THE 4 TRANSCRIPTOMES AND WRITING THE COMMAND-LINE USED 


# Added this so SPAdes knows we're still in our directory and that it can create the assembly files here
current_directory = os.path.dirname(__file__)

# To change the current working directory to the directory of the script
os.chdir(current_directory)

# To define the output directory for SPAdes'assembly
spades_output_dir = os.path.join(current_directory, "ASSEMBLY")

# The SPAdes command based on our class PowerPoint *could also be run in the terminal if desired
spades_cmd = f"spades.py -o ASSEMBLY --s1 30sample_reads_1_bowtie2.fastq --s2 33sample_reads_1_bowtie2.fastq --s3 44sample_reads_1_bowtie2.fastq --s4 45sample_reads_1_bowtie2.fastq"

try:
    subprocess.run(spades_cmd, shell=True, check=True)
    print("SPAdes assembly completed successfully.")
except subprocess.CalledProcessError as e:
    print(f"Error running SPAdes: {e}")

def write_spades_command_to_log(command, log_file):
    with open(log_file, "a") as log:  # Will open the log file in append mode
        log.write("SPAdes command:\n")
        log.write(command + "\n")

# To define the SPAdes command
spades_cmd = "spades.py -o ASSEMBLY --s1 30sample_reads_1_bowtie2.fastq --s2 33sample_reads_1_bowtie2.fastq --s3 44sample_reads_1_bowtie2.fastq --s4 45sample_reads_1_bowtie2.fastq"

# To write the command to the log file
write_spades_command_to_log(spades_cmd, "PipelineProject.log")


# To solve PART 4, CALCULATING NUMBER OF CONTIGS
def calculate_contig_stats(assembly_file): # To calculate the number of contigs that are longer than 1000 base pairs and the total of base pairs contained within the contigs from the assembly file
    num_contigs_over_1000bp = 0 # To initialize the variables to 0. Ready for counting!
    total_bp_over_1000bp = 0

    with open(assembly_file, 'r') as assembly: # To open the assembly file in read mode
        for line in assembly: # To iterate over each line in the assembly 
            if line.startswith('>'): # To identify the header of a fasta file which starts with a >
                header = line.strip().split('_') # To split whitespace
                try:
                    length_index = header.index('length') # This will look for the index of the substring to figure out the length of the contig
                    length = int(header[length_index + 1]) # The formula of the above
                    if length > 1000: # If the contig length if over 1000 bp it will count it 
                        num_contigs_over_1000bp += 1  
                        total_bp_over_1000bp += length # The count will be added here
                except ValueError:
                    print("Error: Unable to extract length from header line:", line.strip())
                except IndexError:
                    print("Error: Index out of range in header line:", line.strip())

    return num_contigs_over_1000bp, total_bp_over_1000bp

def retrieve_longest_contig(assembly_file, output_file): # This will identify the longest contig from the assembly file and save it to a different file
    longest_contig_name = "" # To store the header line of the longest contig found
    longest_contig_length = 0 # To store the length of it

    with open(assembly_file, 'r') as assembly: # To open the file I name later as ASSEMBLY where the assemby will be written
        for line in assembly: # The rest chunk of code is similar to before when counting contigs
            if line.startswith('>'):
                header = line.strip().split('_')
                try:
                    length_index = header.index('length')
                    length = int(header[length_index + 1])
                    if length > longest_contig_length:
                        longest_contig_length = length
                        longest_contig_name = line.strip()
                except ValueError:
                    print("Error: Unable to extract length from header line:", line.strip())
                except IndexError:
                    print("Error: Index out of range in header line:", line.strip())

    with open(assembly_file, 'r') as assembly:
        with open(output_file, 'w') as output:
            for line in assembly:
                if line.strip() == longest_contig_name:
                    output.write(line)
                    for line in assembly:
                        if line.startswith('>'):
                            break
                        output.write(line)
                    break
# To solve PART 5, RETRIEVE LONGEST CONTIG, BLAST+ AGAINST HERPES SUBFAMILY DATABASE, OUTPUT BEST ALIGNMENT 10 HITS WITH HEADER ROW
def perform_blast():
    query_file = "LONGEST_CONTIG" # This will define the retrieved longest contig from the SPAdes assembly
    database_name = "betaherpesvirinae" # This will define the database so we can Blast against it *More instructions on how to do it in my repo
    log_file = "PipelineProject.log" # Tho define where all answers will be written to

    blast_cmd = [ # To define the Blast+ command
        "blastn",  # Chose blastn since nucleotide-nucleotide are being compared
        "-query", query_file,  # Query already defined as the longest contig
        "-db", database_name,  # Database already defined as betaherpesvirinae (subfamily)
        "-outfmt", "6 sacc pident length qstart qend sstart send bitscore evalue stitle",  # Output format
        "-max_target_seqs", "10",  # The number of hits wanted
        "-max_hsps", "1"  # To keep only the best alignment, high scoring pairs
    ]

    try:
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True) # To run BLAST+

        with open(log_file, "a") as f:
            f.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n") # To write the header row as wanted
            f.write(result.stdout) # To write the top hits only
    except subprocess.CalledProcessError as e:
        print(f"Error running Blast+: {e}")
        return

# The precious Main Function 
if __name__ == "__main__":
    if len(sys.argv) < 10:
        print("Usage: python script.py <index> <read1_1> <read1_2> <read2_1> <read2_2> <read3_1> <read3_2> <read4_1> <read4_2>") # The command-line expected to run this python wrapper
        sys.exit(1)

    index_prefix = sys.argv[1] # So the code assumes that the first word in the command-line is the script itself
    read_files = sys.argv[2:] # This will assume the rest are the index and then sample files to use
    
    read_pairs_before, read_pairs_after = process_paired_end_reads(read_files, index_prefix) # To process the read files with the actual Bowtie2 filtering
    
    log_filename = "PipelineProject.log" # The .log file that will contain all of the answers
    
    write_to_log(log_filename, read_pairs_before, read_pairs_after) # For the Bowtie2 part

    spades_output_dir = "ASSEMBLY" # For the SPAdes part

    spades_cmd = f"spades.py -o {spades_output_dir} --s1 30sample_reads_1_bowtie2.fastq --s2 33sample_reads_1_bowtie2.fastq --s3 44sample_reads_1_bowtie2.fastq --s4 45sample_reads_1_bowtie2.fastq"

    write_spades_command_to_log(spades_cmd, log_filename) 

    try:
        subprocess.run(spades_cmd, shell=True, check=True)
        print("SPAdes assembly completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error running SPAdes: {e}")
    num_contigs_over_1000bp, total_bp_over_1000bp = calculate_contig_stats(os.path.join(spades_output_dir, "contigs.fasta")) # To calculate the contigs
    
    with open(log_filename, "a") as log_file: # To write the contig statistics to the log file
        log_file.write(f"There are {num_contigs_over_1000bp} contigs > 1000 bp in the assembly.\n")
        log_file.write(f"There are {total_bp_over_1000bp} bp in the assembly.\n")
    
    longest_contig_file = "LONGEST_CONTIG" # To retrieve and save the longest contig
    retrieve_longest_contig(os.path.join(spades_output_dir, "contigs.fasta"), longest_contig_file)

    perform_blast() # To BLAST

    # Just to print messages of success...
    print(f"Bowtie2 counts before/after filtering written to {log_filename}.")
    print(f"SPAdes assembly completed. Command-line written to {log_filename}.")
    print(f"Assembly contigs written to {log_filename}")
    print(f"Longest contig saved in {longest_contig_file}.")
    print(f"Pipeline complete! See {log_filename}.")

# Done! This code was tested many times with the files provided in this repo. Run time is ~ 2 ish minutes.
# Answers in PipelineProject.log file should look something like this,
# Donor 1 had 40000 read pairs before Bowtie2 filtering and 27103 read pairs after.
# Donor 2 had 40000 read pairs before Bowtie2 filtering and 23569 read pairs after.
# SPAdes command:
# spades.py -o ASSEMBLY --s1 30sample_reads_1_bowtie2.fastq --s2 33sample_reads_1_bowtie2.fastq --s3 44sample_reads_1_bowtie2.fastq --s4 45sample_reads_1_bowtie2.fastq
# There are 44 contigs > 1000 bp in the assembly.
# There are 194635 bp in the assembly.
# sacc	pident	length	qstart	qend	sstart	send	bitscore	evalue	stitle
# NC_006273.2	99.621	7924	9944	17865	36133	28214	14462	0.0	NC_006273.2 Human herpesvirus 5 strain Merlin, complete genome
# NC_003521.1	79.863	1465	11150	12594	32592	31148	1035	0.0	NC_003521.1 Panine herpesvirus 2 strain Heberling, complete genome






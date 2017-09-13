##--------------------------------------------------------------------------------
##
## Final Project: Comparing rRNA depletion methods
##
## Class: PCE Data Science Methods Class
##
## Name: Ben Larson
##
##--------------------------------------------------------------------------------

# Change working directory to project directory
# setwd("/data/basespace")

## Load Libraries----
library(logging)
library(dplyr)
library(stringr)
library(Rsamtools)

## Declare functions

# Set up a logger that goes to a file but suppresses console output
createLog <- function(logFile, loggingLevel) {
    # Remove the test file (invisibly)
    suppressWarnings(invisible(file.remove(logFile)))
    
    # Reset the logging
    logReset()
    
    # Create the logger
    x <- getLogger()
    
    # Set the level of the logger
    setLevel(loggingLevel, x)
    
    # Create a handler, set it's level, and direct it to the logger
    # NOTE: Change first argument to writeToConsole to show on screen
    addHandler(writeToFile, file = logFile, level = loggingLevel, logger = x)
    
    # Return the logger
    return(x)
}

# Function to index the reference genome store in a .fasta file
generate_index <- function(ref.genome) {
    # Construct STAR command with system-specific attributes
    generate.index <- paste0(star.binary,
                             " --runMode genomeGenerate",
                             " --genomeDir ", genome.dir,
                             " --genomeFastaFiles ", ref.genome)
    
    # Execute commnd in the system environment
    system(generate.index)
    loginfo(paste("Generated genome index based on", ref.genome))
}

# Function to align the sequence data with STAR
align_data <- function(data.file1, data.file2, output.prefix,
                       num.thread, read.cmd) {
    # Construct the output path and file prefix
    path.plus.prefix <- paste0(star.output, output.prefix)
    
    # Construct STAR command
    # NOTE: last line eliminates the need to run the samtools "sort" command.
    # --> True, but memory requirements are ridiculous, use samtools instead
    align.data <- paste0(star.binary,
                         " --runThreadN ", num.thread,
                         " --genomeDir ", genome.dir,
                         read.cmd,
                         " --readFilesIn ", data.file1, " ", data.file2, 
                         " --outFileNamePrefix ", path.plus.prefix,
                         " --outSAMtype BAM Unsorted")
   
    # Execute commnd in the system environment
    loginfo(paste("Started mapping", output.prefix))
    system(align.data)
    loginfo(paste("Finished mapping", output.prefix))
    
    # Construct the samtools sort commmand
    sort.data <- paste0("samtools sort",
                        " -m 4G",
                        " -o ", path.plus.prefix, "Aligned.sortedByCoord.out.bam",
                        " -@ ", num.thread, 
                        " ", path.plus.prefix, "Aligned.out.bam")
    
    # Execute commnd in the system environment
    loginfo(paste("Started sorting", output.prefix))
    system(sort.data)
    loginfo(paste("Finished sorting", output.prefix))
}

# Function to generate a pile-up data frame from the results of the mapping 
# NOTE: This is basically a histogram of mapped reads for all positions in the
#       reference genome
make_pile_up <- function(sorted.bam) {
    # Index the .bam file with Rsamtools to create .bai file
    # Returns the name of the .bai file created in the file system
    loginfo(paste(".bai file name: ",indexBam(sorted.bam)))
    
    # Create pileup data frame from .bai file
    # NOTE: The pileup command takes in the sorted .bam file name and looks for a 
    #       .bai file with the same name (including .bam extesion) as the sorted .bam
    #       file, so .bai extension does NOT need to be specified
    x <- pileup(sorted.bam)
    
    # Group total counts by position 
    x <- group_by(x, pos) %>% summarise(tot.count = sum(count))
    return(x)
}

## Main Script----
if(interactive()) {
    # Start the logger
    log.me <- createLog("basespace.remote.log", "DEBUG")
    
    # Specify project directory where all project files and sub directories are
    # stored.
    project.dir <- "/data/basespace/"
    
    # Specify genome subdirectory
    genome.dir <- paste0(project.dir, "genome_index/")
    
    # Specify data directory containing sample data
    data.dir <- paste0(project.dir,"fastq_data/")
    
    # Specify output directory where mapped files will be stored
    star.output <- paste0(project.dir,"star_output/")
    
    # Specify location of STAR binary for mapping sequence data to genome
    # star.binary <- "/home/ec2-user/STAR/bin/Linux_x86_64_static/STAR"
    star.binary <- "~/STAR/source/STAR"
    
    # Specify the fasta file for refence genome
    fasta.file <- paste0(project.dir, "rRNA.fasta")
    
    # Generate genome index files
    # NOTE: only needed once for a given reference
    generate_index(fasta.file)
    
    # Specify paired files with sequence data
    sample.name1 <- list("L11_S1_L001_R1_001.fastq.gz",
                         "L11_S1_L002_R1_001.fastq.gz",
                         "L12_S1_L001_R1_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L001_R1_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L003_R1_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L004_R1_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L001_R1_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L003_R1_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L004_R1_001.fastq.gz",
                         "L12_S1_L002_R1_001.fastq.gz",
                         "L14_S2_L001_R1_001.fastq.gz",
                         "L14_S2_L002_R1_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L002_R1_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L002_R1_001.fastq.gz")
    
    sample.name2 <- list("L11_S1_L001_R2_001.fastq.gz",
                         "L11_S1_L002_R2_001.fastq.gz",
                         "L12_S1_L001_R2_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L001_R2_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L003_R2_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L004_R2_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L001_R2_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L003_R2_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L004_R2_001.fastq.gz",
                         "L12_S1_L002_R2_001.fastq.gz",
                         "L14_S2_L001_R2_001.fastq.gz",
                         "L14_S2_L002_R2_001.fastq.gz",
                         "RED07-10ng-Depleted_S2_L002_R2_001.fastq.gz",
                         "RED08-200ng-Depleted_S1_L002_R2_001.fastq.gz")
    
    # Specify output file name prefix for each pair of files
    output.prefix <- list("L11_L1_",
                          "L11_L2_",
                          "L12_L1_",
                          "Red7_L1_",
                          "Red7_L3_",
                          "Red7_L4_",
                          "Red8_L1_",
                          "Red8_L3_",
                          "Red8_L4_",
                          "L12_L2_",
                          "L14_L1_",
                          "L14_L2_",
                          "Red7_L2_",
                          "Red8_L2_")
    
    # Convert name to full file path
    data.file1 <- lapply(sample.name1, function(x) {paste0(data.dir,x)})
    data.file2 <- lapply(sample.name2, function(x) {paste0(data.dir,x)})
    
    # Set file read command
    # If .fastq files are compressed
    read.cmd <- " --readFilesCommand zcat"
    # If .fastq files are uncompressed
    # read.cmd <- ""
    
    # Call function to align sequence data
    # NOTE: result is a sorted .bam file in the file system
    mapply(FUN = align_data, data.file1, data.file2, output.prefix,
           MoreArgs = list(num.thread = 62, read.cmd = read.cmd))
    
    # Construct sorted .bam file name
    sorted.bam <- lapply(output.prefix, FUN = 
                         function(x) {paste0(star.output, x,
                                             "Aligned.sortedByCoord.out.bam")})
    
    # Call function to create the pile-up
    pup <- lapply(sorted.bam, make_pile_up)
    # Add names to make results distinguishable
    names(pup) <- output.prefix
    
    # Plot a pile-up
    # plot(pup$pos, pup$tot.count)
    
    # No filename specified, default = ".Rdata"
    save(pup, file = "larson.Rdata")
    loginfo("Saved file as larson.Rdata")
    
    logdebug("Make sure you download the data file from the remote host")
    logdebug("IP Address as of June 8, 2017: <ip address>")
    logdebug("Following scip command assumes key file (.pem) is in ~/.ssh/")
    logdebug(paste0("scp -i ~/.ssh/MyKeyPair.pem ec2-user@<ip address>:",
                    "/data/basespace/larson.Rdata ",
                    "<your directory>"))
    
    print(paste("***ATTENTION***"))
    print(paste("You ***MUST*** Execute the scp command below with",
                "<your directory> to run the rest of the code"))
    print(paste0("scp -i ~/.ssh/MyKeyPair.pem ec2-user@<ip address>:",
                 "/data/basespace/larson.Rdata ",
                 "<your directory>"))
}
#---------------------------------------------------------------------------------
# --> Stop remote work and download files to work locally
#---------------------------------------------------------------------------------
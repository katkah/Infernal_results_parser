#!/bin/bash

#PBS -l select=1:ncpus=6:mem=35gb:scratch_local=100gb
#PBS -l walltime=20:00:00
#PBS -N download
#

# Script uses edirect utilities to download assembly genomes based on the query it outputs taxonomy.csv

######################################################################################
#SET VARIABLES IN THIS SECTION

# Set the output CSV file
output_csv="taxonomy.csv"

# Set the query
#query='"Oxalidales"[Organism] AND (latest[filter] AND all[filter] NOT anomalous[filter])'
query='"Coleoptera"[Organism] AND "representative genome"[filter]'


# Set the output directory
output_directory="full/path/to/database/"
######################################################################################

# Create output_directory
mkdir -p $output_directory
cd $output_directory

# Load edirect utilities
module add edirect

# Prepare two temporary file with taxid and ftp links by searching ncbi assembly database
esearch -db assembly -query "${query}" | esummary | xtract -pattern DocumentSummary -sep ';' -element Taxid,FtpPath_GenBank > taxid_line.txt

# Taxonomy ranks to extract
taxonomy_ranks=("Superkingdom" "Kingdom" "Phylum" "Class" "Order" "Family" "Genus")

# Add the header to CSV
echo "superkingdom,kingdom,phylum,class,order,family,genus,ScientificName,GCA,taxid" > "${output_csv}"

#Notes:
#Examples of how to use edirect utilities:
#this command will give you the ScientificName for id 138855
#efetch -db taxonomy -id 138855 -format xml | xtract -pattern Taxon -element ScientificName
#this command fetches a selected rank (- "kingdom" in this case) for id 13885
#esearch -db taxonomy -query "138855 [taxID]" | efetch -format native -mode xml | xtract -pattern Taxon -block "*/Taxon" -if Rank -equals "kingdom" -tab "\n" -element ScientificName


# Step 1: Download assemblies in FASTA format and retrieve taxonomy information

for row in $(cat taxid_line.txt); do
    taxid=$(echo "$row" | cut -d ';' -f 1)
    line=$(echo "${row#*;}")
    
    echo "Processing line: $line"
    echo "TaxID: $taxid"
    
    # It greps only GCA file name from line e.g. ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/024/708/375/GCA_024708375.1_oxalis_stricta_v0.1
    fname=$(echo $line | grep -o 'GCA_.*' | sed 's/$/_genomic.fna.gz/') ;

    echo "Assembly file: $fname"

    # Fetch taxonomy_ranks information using Taxid
    echo "Fetching taxonomy info:"
    for rank in "${taxonomy_ranks[@]}"; do
        taxonomy_info=$(esearch -db taxonomy -query "${taxid} [taxID]" | efetch -format native -mode xml | xtract -pattern Taxon -block "*/Taxon" -if Rank -equals "${rank}" -element ScientificName)

        # Check if taxonomy_info is empty and set it to "NA" if true; that means the download was not successful even after the x-th attempt
        max_attempts=100
        attempts=0

        while [ "$attempts" -lt "$max_attempts" ]; do
            # If taxonomy_info length not zero then fetching info was success and break the loop
            if [ -n "$taxonomy_info" ]; then
                break
            fi
            taxonomy_info=$(esearch -db taxonomy -query "${taxid} [taxID]" | efetch -format native -mode xml | xtract -pattern Taxon -block "*/Taxon" -if Rank -equals "${rank}" -element ScientificName)
            attempts=$((attempts + 1))
        done
        # If taxonomy_info length is zero then it was unsuccess even after the x-th attempt and set it to NA
        if [ -z "$taxonomy_info" ]; then
            taxonomy_info="NA"
            echo "fetching taxonomy_info for rank ${rank} and ${taxid} was unsuccessful"
        fi

        # Echo without the newline character at the end
        echo -n "${taxonomy_info}," >> "${output_csv}"
    done
    
    # Fetch the "ScientificName" information using Taxid
    taxon=$(efetch -db taxonomy -id "${taxid}" -format xml | xtract -pattern Taxon -element ScientificName)
    max_attempts=100
    attempts=0

    while [ "$attempts" -lt "$max_attempts" ]; do
        if [ -n "$taxon" ]; then
            break
        fi
        taxon=$(efetch -db taxonomy -id "${taxid}" -format xml | xtract -pattern Taxon -element ScientificName)
        attempts=$((attempts + 1))
    done    
    
    
    # Check if taxon is empty and set it to "NA" if true       
    if [ -z "$taxon" ]; then
        taxon="NA"
        echo "fetching taxon for and ${taxid} was unsuccessful"
    fi
    # Echo without the newline character at the end
    echo -n "${taxon}," >> "${output_csv}"
    
    # Echo GCA with the newline character at the end
    GCA=$(echo "${fname}" | cut -d '_' -f2)
    echo -n "GCA_${GCA}," >> "${output_csv}"
    echo "$taxid" >> "${output_csv}"
    
    # Downloading assembly file
    max_attempts=3
    attempts=0

    while [ "$attempts" -lt "$max_attempts" ]; do
        wget "$line/$fname"

        if [ $? -eq 0 ]; then
            echo "File downloaded successfully"
            break  # Break out of the loop if download is successful
        else
            echo "Error: File download failed, retrying..."
            attempts=$((attempts + 1))
        fi
    done
    
    if [ "$attempts" -ge "$max_attempts" ]; then
        echo "Error: Maximum download attempts reached. File download unsuccessful."
    fi
    echo "Finished processing line."
done

# Remove temporary files
rm taxid_line.txt


# Step2: Unzip downloaded files
for file in *genomic.fna.gz
do
    gunzip "$file"
done



echo "Download and CSV creation complete. Check the '${output_directory}' directory for downloaded assemblies and the '${output_csv}' file for information."


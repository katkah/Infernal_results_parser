import os
import csv
import pandas as pd
import sys
import subprocess

def reverse_complement(dna_string):
    complement_dict = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}

    reversed_dna = dna_string[::-1]
    reverse_complement_dna = ''

    for base in reversed_dna:
        if base.upper() in complement_dict:
            complement_base = complement_dict[base.upper()]
            reverse_complement_dna += complement_base if base.isupper() else complement_base.lower()
        else:
            reverse_complement_dna += base

    return reverse_complement_dna


def process_genomic_files(database_dir, result_dir):
    # Get a list of genomic.bed files in result_dir
    bed_files = [file for file in os.listdir(result_dir) if file.endswith("genomic.bed")]

    # Dictionary to store the results
    results = {}
    
    # Iterate over each genomic.bed file
    for bed_file in bed_files:
        # Construct the corresponding genomic.fna file name
        genomic_fna_file = bed_file.replace("genomic.bed", "genomic.fna")
        

        
        # Build the full paths for input and output files
        input_genomic_fna = os.path.join(database_dir, genomic_fna_file)
        #output_genomic_fasta = os.path.join(result_dir, bed_file.replace("genomic.bed", "genomic.fasta"))

        # Construct the subprocess command
        
        cmd = f"seqtk subseq {input_genomic_fna} {os.path.join(result_dir, bed_file)}"
        
        print(f"running {cmd}")
        # Execute the subprocess command and capture the output
        process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = process.communicate()

        # Check if the command was successful
        if process.returncode == 0:
            # Store the output in the results dictionary
            results[bed_file] = output.decode('utf-8')
        else:
            # Print an error message if the command failed
            print(f"Error processing {bed_file}: {error.decode('utf-8')}")
    
    return results
"""
def parse_seqtk_results(results, df, extend_region):
    for k in results.keys():
        full_gca_name = k.replace(".bed", ".fna")
        records = results[k].split('>')
        # If the file starts with '>', the first record will be empty
        if records[0] == '':
            records = records[1:]
        # There is a header and a sequence in each record. The sequence might be spread on more lines divided by "\n". The sequence might be also missing.
        for record in records:
            parts = record.strip().split("\n")
            header = parts[0]
            #sequence parts[1:]
            if len(p) == 1:
                seq = "Unable to extract sequence"
            else:
                seq = ''.join(parts[1:])                        
            target_name = header.split(":")[0]
            coordinates = header.split(":")[1]
            seq_from_extend = coordinates.split('-')[0]
            seq_to_extend = coordinates.split('-')[1]
            
            #TODO make columns seq_from_extend, seq_to_extend before 

"""


def open_and_extract_taxonomy(file_path):
    """
    Open and extract taxonomy file.

    Args:
        file_path (str): The path to the taxonomy file.

    Returns:
        DataFrame: A DataFrame containing taxonomy.
    """
    # Check if the file exists
    if os.path.exists(file_path):
        # File exists, so read CSV file into a DataFrame
        df = pd.read_csv(file_path)
        # Drop the first column. The first column is just numbers => it's an artefact from R script 02
        df = df.drop(df.columns[0], axis=1)

        return df
    else:
        print(f"The file '{file_path}' does not exist.")
        return

def parse_genomic_file_name(file_name):
    """
    Parse the name of a genomic file based on underscores. e.g result_U1_vs_GCA_008126665.1_ASM812666v1_genomic.csv

    Args:
        file_name (str): The name of the genomic file.

    Returns:
        dict: A dictionary containing parsed components of the file name.
    """
    # Split the file name based on underscores
    components = file_name.split('_')
    full_name = components[3:]
    full_name = '_'.join(full_name)
    full_name = full_name.replace("genomic.csv", "genomic.fna")

    # Extract relevant information based on the position of elements
    result = {
        'model': components[1],
        'assembly': "GCA_" + components[4],
        'full_GCA_name': full_name
    }

    return result

def extract_info_from_genomic_file(file_path,parsed_info):
    """
    Parse the content of a genomic file based on white spaces.

    Args:
        file_path (str): The name of the genomic file. The genomic file is not tab separated, there is unspecified number of spaces between columns.
        parsed_info (dict): Dictionary of two fields extracted from the file name e.g. {'model': 'U1', 'assembly': 'GCA_015342455.1'}

    Returns:
        DataFrame: A DataFrame containing all infernal matches found in the file.
    """
    
    print(f"Processing file: {file_path}")
    with open(file_path, 'r') as file:
        df_header = ["target_name","accession1","query_name","accession2","mdl","mdl_from","mdl_to","seq_from","seq_to","strand","trunc","pass","gc","bias","score","E-value","inc"]
        
        # Create an empty DataFrame to store the data
        df = pd.DataFrame(columns=df_header)                       
        
        # Read the file line by line, skipping lines starting with # - those are comments
        for line in file:
            if not line.startswith('#'):
                #Take first 17 fields and add them to the DataFrame. Omit the last column - "description of target". The column is problematic, because it contains spaces. Fortunately, it is not important and can be omitted. 
                important_fields = line.split()[:17]
                df.loc[len(df.index)] = important_fields


    #Add columns 'model' and 'assembly' to the DataFrame
    df['model'] = parsed_info['model']
    df['GCA'] = parsed_info ['assembly']
    df['full_GCA_name'] = parsed_info['full_GCA_name']
    df['number'] = pd.NA

    return df   

def create_string(row):
    return f"{row['target_name']}/{row['seq_from']}-{row['seq_to']}"
    
def extract_info_from_alignment_file(alignment_file, df):
    print(f"Processing file: {alignment_file}")
    # Extract alignment for each hit; The line of interest starts with {target_name}/{seq_from}-{seq_to}; The "Maybe hit" designated with a "?" are not listed in alignment files. 
    for index, row in df.iterrows():
        search_string = create_string(row)
        with open(alignment_file, 'r') as file:
            for line in file:              
                if line.startswith(search_string):
                # Split the line based on search_string and save the second part as 'alignment'
                    alignment_value = line.split(search_string)[1].strip()
                    df.at[index, 'alignment'] = alignment_value
                    
    
    # Extract consensus structure and consensus alignment   
    search_string = "#=GC SS_cons"
    search_string2 = "#=GC RF"
    structure = pd.NA
    gc_rf = pd.NA

    with open(alignment_file,'r') as file:
        for line in file:              
            if line.startswith(search_string):
            # Split the line based on search_string and save the second part as 'SS_cons'
                structure = line.split(search_string)[1].strip()

            # Split the line based on search_string2 and save the second part as 'GC_RF'
            if line.startswith(search_string2):
                gc_rf = line.split(search_string2)[1].strip()

    
    df['GC_SS_cons'] = structure
    df['GC_RF'] = gc_rf
    
    return df           
        

        
def extract_genomic_ali_files(directory_path):
    """
    Parse the content of all genomic.csv and genomic-alignment files.

    Args:
        directory_path (str): The name of the directory containing genomic files to be processed.

    Returns:
        DataFrame: A DataFrame containing concatenated infernal matches found in all the genomic files. Plus structure and alignment from genomic-alignment files. 
    """    
    
    # Check if the directory exists
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' does not exist.")
        sys.exit(1)
    
    # List all files in the directory
    files = os.listdir(directory_path)
    
    # Filter files that end with "_genomic.csv"
    genomic_files = [file for file in files if file.endswith('_genomic.csv')]
    alignment_files = [file for file in files if file.endswith('genomic-alignment')]
    
    
    if not genomic_files:
        print(f"No files with the '_genomic.csv' extension found in '{directory_path}'.")
        return
    if not alignment_files:
        print(f"No files with the 'genomic-alignment' extension found in '{directory_path}'.")
        return

    # Open and extract information from each genomic file and correspondent alignment file, save it to two dataframes (hits and maybe)
    df_hits_all = None
    df_maybe_all = None
    #Save number of hits and maybe for each GCA file to report at the end
    gca_number_hits = {'full_GCA_name': [], 'number_of_hits' : [], 'number_of_possible_hits':[]} 
 
    for file_name in genomic_files:
    
        # Extract info from genomic file name e.g. {'model': 'U1', 'assembly': 'GCA_015342455.1'}
        parsed_info = parse_genomic_file_name(file_name)
        
        #Extract content of a genomic file
        file_path = os.path.join(directory_path, file_name)
        df_file = extract_info_from_genomic_file(file_path, parsed_info)
        
        # Divide the DataFrame into two based on the "inc" column
        df_hits = df_file[df_file['inc'] == '!'].copy()
        df_maybe = df_file[df_file['inc'] != '!'].copy()
      
        #Extract content of alignment file for df_hits     
        alignment_file = file_name.replace("genomic.csv", "genomic-alignment")
        alignment_path = os.path.join(directory_path, alignment_file)
        if alignment_file not in alignment_files:
            print(f"file {alignment_file} does not exist!")
        else:
            df_hits = extract_info_from_alignment_file(alignment_path, df_hits)
        
        # Count number of hits
        gca_number_hits['full_GCA_name'].append(parsed_info['full_GCA_name'])
        gca_number_hits['number_of_hits'].append(len(df_hits))       
        gca_number_hits['number_of_possible_hits'].append(len(df_maybe))
     
        # Sorting the DataFrame by the "E-value" column in descending order
        df_hits.sort_values(by='E-value', ascending=False, inplace=True)    
        df_maybe.sort_values(by='E-value', ascending=False, inplace=True) 
               
        # Change "number" with order information
        for index,  row in df_hits.iterrows():
            df_hits.at[index, 'number'] = index + 1
        
        # Number can't be int because of NaN values in the column      
        for index, row in df_maybe.iterrows():
            df_maybe.at[index, 'number'] = index + 1 
        
        # Concatenate df_hits_all and df_hits along rows
        if df_hits_all is not None:
            df_hits_all = pd.concat([df_hits_all, df_hits], axis=0)
        # If the first file is being parsed, create df 
        else:
            df_hits_all = df_hits 

        # Concatenate df_maybe_all and df_maybe along rows
        if df_maybe_all is not None:
            df_maybe_all = pd.concat([df_maybe_all, df_maybe], axis=0)
        # If the first file is being parsed, create df 
        else:
            df_maybe_all = df_maybe
        
    gca_number_hits = pd.DataFrame(gca_number_hits)    
    return df_hits_all, df_maybe_all, gca_number_hits



"""
def add_extended_region(df,directory_database_path):
    for index, row in df.iterrows():
        
        extended_region = pd.NA
        
        target = f">{row['target_name']}"
        
        assembly_file = row['full_GCA_name']
        assembly_file = os.path.join(directory_database_path, assembly_file)
        
        sequence = extract_sequence(assembly_file, target)
        
        if row['strand'] == "+":
            extended_region = sequence[row['new_start']:row['new_end']]
        if row['strand'] == "-":
            extended_region = sequence[row['new_end']:row['new_start']]
        
        df.at[index, 'extended_region'] = extended_region     
    
    return df
"""
"""
#Function creates coordinates.bed one GCA file, saves it to result_dir
def save_dataframe_as_bed(bed_filename, df, result_dir, extend_region):


    bed_filepath = os.path.join(result_dir, bed_filename)
    
    # Create a copy of the selected columns to avoid modifying the original DataFrame
    selected_df = df.copy()

    #Convert columns to numeric type
    selected_df['seq_from'] = pd.to_numeric(selected_df['seq_from'], errors='coerce')
    selected_df['seq_to'] = pd.to_numeric(selected_df['seq_to'], errors='coerce')
    
    # Check the strand column and swap seq_from and seq_to if strand is minus
    minus_strand_mask = selected_df['strand'] == '-'
    selected_df.loc[minus_strand_mask, ['seq_from', 'seq_to']] = selected_df.loc[minus_strand_mask, ['seq_to', 'seq_from']].values
    
    # Subtract 1 from the 'seq_from' column for seqtk to work properly
    selected_df['seq_from'] = selected_df['seq_from'] - 1 - extend_region
    if selected_df['seq_from'] < 0:
        selected_df['seq_from'] = 0
    selected_df['seq_to'] = selected_df['seq_to'] + extend_region
       
    
    # Save the dataframe to the bed file without header and columns tab delimited
    selected_df.to_csv(bed_filepath, header=False, index=False, sep='\t')
    

    
    return 
"""   
"""
#Function creates coordinates.bed for all GCA files
def process_gca_files(df, directory_database_path, result_dir, extend_region ):
    # Create list of unique GCA files    
    GCA_list =  df['full_GCA_name'].unique()    
    # Get targets for each GCA file
    for gca in GCA_list:        

        # Check if directory_database_path exists
        if not os.path.exists(directory_database_path):
            print(f"Directory '{directory_database_path}' does not exist.")
            sys.exit(1)
        
        # Check if GCA file exists
        files = os.listdir(directory_database_path)       
        gca_files = [file for file in files if file.endswith('genomic.fna')]   
        file_path = os.path.join(directory_database_path, gca)      
        if gca not in gca_files:
            print(f"file {gca_file} does not exist!")
        else:
            print(f"Processing file: {file_path}")
            selected_columns = ["target_name","seq_from","seq_to", "strand", "full_GCA_name"]            
            targets = df.loc[df['full_GCA_name'] == gca,selected_columns]
            full = df.loc[df['full_GCA_name'] == gca, ["targed_name"]]
            # Change "genomic.fna" to "genomic.bed" in the full_GCA_name
            
            save_dataframe_as_bed(gca.replace("genomic.fna", "genomic.bed"), targets, result_dir, extend_region)
            save_dataframe_as_bed(gca.replace(".fna", "_full.bed"), targets, result_dir, 0) #save whole sequence to count the length
                
    return     
"""

def read_result(seqtk_result):
    parts = seqtk_result.strip().split("\n")
    header = parts[0]

    if len(parts) == 1:
        seq = "Unable to extract sequence"
    else:
        seq = ''.join(parts[1:])    
                            
    #target_name = header.split(":")[0]
    coordinates = header.split(":")[1]
    seq_from_extend = coordinates.split('-')[0]
    seq_to_extend = coordinates.split('-')[1]
    
    return seq_from_extend, seq_to_extend, seq

def call_seqtk(row, database_dir, result_dir, extend_region):
    gca =  row['full_GCA_name']    
    gca_path = os.path.join(database_dir, gca)
    files = os.listdir(database_dir)
    if gca not in files:
        print(f"file {gca} does not exist!")
        
    bed_filename = "tmp.bed"
    bed_filepath = os.path.join(result_dir, bed_filename)
    
    if row['strand'] == '-':
        seq_from = int(row['seq_to'])
        seq_to = int(row['seq_from'])
    else:
        seq_from = int(row['seq_from'])
        seq_to = int(row['seq_to'])
        
    seq_from = seq_from -1 - extend_region
    if seq_from < 0:
        seq_from = 0
    seq_to = seq_to + extend_region


    data = [row['target_name'], seq_from, seq_to]

    with open(bed_filepath, 'w', encoding='UTF8', newline='') as f:
        writer = csv.writer(f,delimiter = '\t')
        writer.writerow(data)
    
    # Build the full paths for input files
    input_genomic_fna = os.path.join(database_dir, gca)
    
    # Construct the subprocess command
    cmd = f"seqtk subseq {input_genomic_fna} {bed_filepath}"
    
    # Execute the subprocess command and capture the output
    process = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    seqtk_result = ""
    # Check if the command was successful
    if process.returncode == 0:
        # Store the output in the results dictionary
        seqtk_result = output.decode('utf-8')
    else:
        # Print an error message if the command failed
        print(f"Error processing {bed_file}: {error.decode('utf-8')}")
    
    sfrom, sto, s = read_result(seqtk_result)
    # Swap sfrom and sto, reverse complement the sequence
    if row['strand'] == '-':
        a = sfrom
        sfrom = sto
        sto = a
        s = reverse_complement(s)
           
    return sfrom, sto, s



def main():

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python parse_tables.py infernal_results_dir taxonomy_file region_extend_length database_sequences_dir path_where_to_save_results")
        sys.exit(1)
        
    # Extracting command line arguments
    directory_path = sys.argv[1]
    taxonomy_file = sys.argv[2]
    extend_region = int(sys.argv[3])
    database_dir = sys.argv[4] 
    result_dir = sys.argv[5]
    result_file_hits = os.path.join(result_dir, "infernal_result_hits.csv")
    result_file_maybe = os.path.join(result_dir, "infernal_result_possible_hits.csv")
    result_file_gca = os.path.join(result_dir, "infernal_result_report.csv")

    if not os.path.exists(database_dir):
        print(f"Directory '{database_dir}' does not exist.")
        sys.exit(1) 

    df_hits, df_maybe, gca_number_hits = extract_genomic_ali_files(directory_path)

    # Merge genomic file with taxonomy based on GCA
    df_taxonomy = open_and_extract_taxonomy(taxonomy_file)
    df_hits = pd.merge(df_hits, df_taxonomy, on='GCA', how='left', suffixes=('_df', '_taxonomy')).fillna('NA')
    df_maybe = pd.merge(df_maybe, df_taxonomy, on='GCA', how='left', suffixes=('_df', '_taxonomy')).fillna('NA')    
        
    # Call seqtk
    for index, row in df_hits.iterrows():   
        seq_from, seq_to, seq = call_seqtk(row,database_dir, result_dir, extend_region)
        df_hits.at[index, 'seq_from_extend'] = seq_from 
        df_hits.at[index, 'seq_to_extend'] = seq_to
        df_hits.at[index, 'extend_sequence'] = seq
        
    # Merge with gca_number_hits
    df_hits = pd.merge(df_hits, gca_number_hits, on='full_GCA_name', how='right', suffixes=('_df', '_gca')).fillna('NA')
    df_maybe = pd.merge(df_maybe, gca_number_hits, on='full_GCA_name', how='right', suffixes=('_df', '_gca')).fillna('NA')
    
    df_h = df_hits.sort_values(['E-value'], ascending=False).groupby('full_GCA_name')
    df_m = df_maybe.sort_values(['E-value'], ascending=False).groupby('full_GCA_name') 
     
    if df_hits.empty: 
        print('Infernal did not find any hits in any of the analyzed genomes!')
    else:
        df_h.to_csv(result_file_hits, index=False, na_rep='NA')
        df_m.to_csv(result_file_maybe, index=False, na_rep='NA')
    
if __name__ == "__main__":
    main()

#This is the expected output:
#csv file infernal_results.csv
#with fields:
#model,GCA,gc,score,evalue,label,number,suma,mdl_from,ID,seq_from,seq_to,align_seq,model_seq,sec_struct,X1,superkingdom,kingdom,phylum,class,order,family,genus,ScientificName


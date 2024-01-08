import os
import csv
import pandas as pd
import sys


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
    full_name = f"GCA_{components[4:]}"
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
                # Split the line based on sarch_string and save the second part as 'alignment'
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

    # Open and extract information from each genomic file and correspondent alignment file, save it to a dataframe
    df=None
    for file_name in genomic_files:
    
        # Extract info from genomic file name e.g. {'model': 'U1', 'assembly': 'GCA_015342455.1'}
        parsed_info = parse_genomic_file_name(file_name)
        
        #Extract content of a genomic file
        file_path = os.path.join(directory_path, file_name)
        df_file = extract_info_from_genomic_file(file_path, parsed_info)
        
        #Extract content of alignment file
        
        alignment_file = file_name.replace("genomic.csv", "genomic-alignment")
        alignment_path = os.path.join(directory_path, alignment_file)
        if alignment_file not in alignment_files:
            print(f"file {alignment_file} does not exist!")
        else:
            df_file = extract_info_from_alignment_file(alignment_path, df_file)
            
        # Concatenate df and df_file along rows
        if df is not None:
            df = pd.concat([df, df_file], axis=0)
        # If the first file is being parsed, create df 
        else:
            df = df_file 
    return df
    
def get_max_length(row,directory_database_path):
    assembly_file = row['full_GCA_name']
    assembly_file = os.path.join(directory_path, alignment_file)
    target = row['target_name']
    target = f">{taget}"
    
    with open(assembly_file,'r') as file:
        for line in lines:
            if line.startswith(target):
    
   
# Function to calculate new_start and new_end based on the strand
def calculate_new_coordinates(row, region_length,directory_database_path):
    if row['strand'] == '+':
        row['new_start'] = int(row['seq_from']) - region_length
        row['new_end'] = int(row['seq_to']) + region_length
    elif row['strand'] == '-':
        row['new_start'] = int(row['seq_from']) + region_length
        row['new_end'] = int(row['seq_to']) - region_length
    
    # Get how long the sequence is - The end of the sequence is max possible coordinate
    max_length = get_max_length(row,directory_database_path)
        
    #Set new_start or new_end to the beginning if it is beyond plausible coordinates
    if row['new_start'] <= 0:
        row['new_start'] = 1
    if row['new_end'] <= 0:
        row['new_end'] = 1
    if row['new_start'] > max_length:
        row['new_start'] = max_length
    if row['new_end'] > max_length:
        row['new_end'] > max_length
        
    return row
    

# Replace 'your_directory_path' with the actual path of the directory you want to process
def main():

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 6:
        print("Usage: python parse_tables.py path_to_directory_with_infernal_results path_to_taxonomy_file region_extend_length path_to_directory_with_database_sequences name_result_file ")
        sys.exit(1)
        
    # Extracting command line arguments
    directory_path = sys.argv[1]
    taxonomy_file = sys.argv[2]
    region_extend_length = int(sys.argv[3])
    directory_database_path = sys.argv[4] 
    result_file = sys.argv[5]

    df = extract_genomic_ali_files(directory_path)

    df_taxonomy = open_and_extract_taxonomy(taxonomy_file)

    # Merge genomic file with taxonomy based on GCA
    df = pd.merge(df, df_taxonomy, on='GCA', how='left', suffixes=('_df', '_taxonomy')).fillna('NA')

    # Add columns with coordinates of extended regions
    df = df.apply(calculate_new_coordinates, axis=1, args=(region_extend_length,directory_database_path,))
    
    #df = add_extended_region(df) 
        
    if df.empty: 
        print('Infernal did not find any hits in any of the analyzed genomes!')
    else:
        df.to_csv(result_file, index=False, na_rep='NA')
    

if __name__ == "__main__":
    main()

#This is the expected output:
#csv file infernal_results.csv
#with fields:
#model,GCA,gc,score,evalue,label,number,suma,mdl_from,ID,seq_from,seq_to,align_seq,model_seq,sec_struct,X1,superkingdom,kingdom,phylum,class,order,family,genus,ScientificName


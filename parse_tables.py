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
    Parse the name of a genomic file based on underscores.

    Args:
        file_name (str): The name of the genomic file.

    Returns:
        dict: A dictionary containing parsed components of the file name.
    """
    # Split the file name based on underscores
    components = file_name.split('_')

    # Extract relevant information based on the position of elements
    result = {
        'model': components[1],
        'assembly': "GCA_" + components[4],
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

    return df   


def open_and_extract_genomic_files(directory_path):
    """
    Parse the content of all genomic files.

    Args:
        directory_path (str): The name of the directory containing genomic files to be processed.

    Returns:
        DataFrame: A DataFrame containing concatenated infernal matches found in all the genomic files.
    """    
  
    # Check if the directory exists
    if not os.path.exists(directory_path):
        print(f"Directory '{directory_path}' does not exist.")
        sys.exit(1)
    
    # List all files in the directory
    files = os.listdir(directory_path)
    
    # Filter files that end with "_genomic.csv"
    genomic_files = [file for file in files if file.endswith('_genomic.csv')]
    
    if not genomic_files:
        print(f"No files with the '_genomic.csv' extension found in '{directory_path}'.")
        return
    
    df=None
    # Open and extract information from each genomic file, save it to a dataframe
    for file_name in genomic_files:
        # Extract info from genomic file name e.g. {'model': 'U1', 'assembly': 'GCA_015342455.1'}
        parsed_info = parse_genomic_file_name(file_name)
        
        file_path = os.path.join(directory_path, file_name)
        df_file = extract_info_from_genomic_file(file_path, parsed_info)
        

        # Concatenating df and df_file along rows
        if df is not None:
            df = pd.concat([df, df_file], axis=0)
        # When the first file is parsed, create df 
        else:
            df = df_file 
    return df

#parsovat argumenty if __mane args
# Replace 'your_directory_path' with the actual path of the directory you want to process
def main():

    # Check if the correct number of arguments is provided
    if len(sys.argv) != 4:
        print("Usage: python parse_tables.py directory_path taxonomy_file result_file")
        sys.exit(1)
        
    # Extracting command line arguments
    directory_path = sys.argv[1]
    taxonomy_file = sys.argv[2]
    result_file = sys.argv[3]       

    df = open_and_extract_genomic_files(directory_path)
    df_taxonomy = open_and_extract_taxonomy(taxonomy_file)

    # Merge genomic file with taxonomy based on GCA
    df = pd.merge(df, df_taxonomy, on='GCA', how='left', suffixes=('_df', '_taxonomy')).fillna(0)

    if df.empty: 
        print('Infernal did not find any hits!')
    else:
        df.to_csv(result_file, index=False)

if __name__ == "__main__":
    main()

#This is the expected output:
#csv file infernal_results.csv
#with fields:
#model,GCA,gc,score,evalue,label,number,suma,mdl_from,ID,seq_from,seq_to,align_seq,model_seq,sec_struct,X1,superkingdom,kingdom,phylum,class,order,family,genus,ScientificName


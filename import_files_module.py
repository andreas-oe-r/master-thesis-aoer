"""
### File name:                                                              
import_files_module.py                              

### General description of the script:                                      
This module contains functions for the import of the main split-vep 
data-file, import of other files containing patient ID, sample ID and 
import of files with gene-specific annotation information.                                                                 
                                                                          
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Imports                                                                     #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# The script uses pandas for storing data in a dataframe object, and to perform
# various basic data analysis steps (extraction, subsetting etc.). 
import pandas as pd

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Functions                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

### Function 0:

def read_patient_ID_file(path_to_patient_ID_file):
    """
    Reads a tab-delimited text file with a single line containing the patient 
    ID.

    Parameters
    ----------
    path_to_patient_ID : str
        The file location of the tab-delimited text file

    Returns
    -------
    patient_id : str
        The patient ID
    """

    # Read file with patient ID
    sample_file = open(path_to_patient_ID_file, "r")
    
    # Read contents of file into a list
    sample_file_list = sample_file.readlines()
    
    # Assign patient ID to variable
    patient_ID = sample_file_list[0]
    
    # Remove trailling newline
    patient_ID = patient_ID.rstrip()
    
    # Close file
    sample_file.close()

    return patient_ID


### Function 1:

def read_project_sample_file(path_to_project_sample_file):
    """
    Reads a tab-delimited text file that contains 2 columns and 90 rows. The 
    first column contains the project sample ID's, and the second column 
    contains the corresponding patient ID.

    Parameters
    ----------
    path_to_project_sample_file : str
        The file location of the tab-delimited text file

    Returns
    -------
    project_sample_ID_list : pandas dataframe object
        A dataframe of dimensions 90 x 2
    """

    project_sample_ID_list = pd.read_csv(
            path_to_project_sample_file,
            sep = '\t',
            header = 0)

    return project_sample_ID_list


### Function 2:

def get_sample_ID_from_patient_ID(patient_ID, project_sample_ID_list):
    """
    Finds the corresponding project_sample_ID for a patient_ID by searching 
    for the patient_ID string in a dataframe containing patient ID's and 
    project sample ID's.

    Parameters
    ----------
    patient_ID : str
        The patient ID
        
    project_sample_ID_list : pandas dataframe object
        A dataframe of dimensions 90 x 2

    Returns
    -------
    project_sample_ID : str
        The project sample ID
    """

    # Convert column type to string  
    project_sample_ID_list['NGC_HPC_patient_ID'] = \
        project_sample_ID_list['NGC_HPC_patient_ID'].astype(str)
        
    # Convert column type to string  
    project_sample_ID_list['Project_sample_ID'] = \
        project_sample_ID_list['Project_sample_ID'].astype(str)

    # Find the row where the patient ID matches the values in the 
    # NGC_HPC_patient_ID column
    project_sample_ID_row = \
        project_sample_ID_list[project_sample_ID_list[
                'NGC_HPC_patient_ID'].str.contains(patient_ID)] 
    
    project_sample_ID = project_sample_ID_row.iloc[0]['Project_sample_ID']

    return project_sample_ID


### Function 3:

def read_cram_file_list(path_to_cram_list_file):
    """
    Reads a tab-delimited text file that contains 2 columns and 90 rows. The 
    first column contains the patient ID, and the second column 
    contains the corresponding path to a cram file.

    Parameters
    ----------
    path_to_cram_list_file : str
        The file location of the tab-delimited text file

    Returns
    -------
    cram_list : pandas dataframe object
        A dataframe of dimensions 90 x 2
    """

    cram_list = pd.read_csv(
            path_to_cram_list_file,
            sep = '\t',
            header = 0)

    return cram_list


### Function 4:

def get_path_to_cram_file_from_patient_ID(patient_ID, cram_file_list):
    """
    Finds the corresponding path to a cram file for a patient_ID by searching 
    for the patient_ID string in a dataframe containing patient ID's and 
    paths to cram files.

    Parameters
    ----------
    patient_ID : str
        The patient ID
        
    cram_file_list : pandas dataframe object
        A dataframe of dimensions 90 x 2

    Returns
    -------
    project_sample_ID : str
        The project sample ID
    """

    # Convert column type to string  
    cram_file_list['Patient_ID'] = \
        cram_file_list['Patient_ID'].astype(str)
        
    # Convert column type to string  
    cram_file_list['path_to_cram_file'] = \
        cram_file_list['path_to_cram_file'].astype(str)

    # Find the row where the patient ID matches the values in the 
    # NGC_HPC_patient_ID column
    crame_file_ID_row = \
        cram_file_list[cram_file_list[
                'Patient_ID'].str.contains(patient_ID)] 
    
    path_to_cram_file = crame_file_ID_row.iloc[0]['path_to_cram_file']

    return path_to_cram_file


### Function 5:

def import_split_vep_file(path_to_split_vep_file_str):
    """
    Reads a tab-delimited text file that contains 100 columns and dynamic 
    number of rows. Rows corresponds to each genomic location where a variant 
    is located. Columns contains information the variant.

    Parameters
    ----------
    path_to_split_vep_file_str : str
        The file location of the tab-delimited text file

    Returns
    -------
    variant_calls_dataframe : pandas dataframe object
        A dataframe of dimensions x x 100
    """
    
    # Import file as a pandas dataframe and add headers for all columns (100)
    variant_calls_dataframe = pd.read_csv(path_to_split_vep_file_str, 
            delim_whitespace = True,
            names = ["Chromosome_Position", "ID", "Reference", "Alternative",
                     "Quality", "Type", "Genotype", "Allele_depth_AD",
                     "Depth_of_coverage_DP", "Phred-scaled_likelihoods_PL",
                     "Genotype_Quality_GQ", "Allele", "Consequence", "Impact",
                     "Symbol", "Gene", "Feature_type", "Feature", "Biotype",
                     "Exon", "Intron", "HGVSc", "HGVSp", "cDNA_position",
                     "CDS_position", "Protein_position", "Amino_acids",
                     "Codons", "Existing_variation", "DISTANCE", "STRAND",
                     "FLAGS", "VARIANT_CLASS", "SYMBOL_SOURCE", "HGNC_ID",
                     "CANONICAL", "MANE", "TSL", "APPRIS", "CCDS", "ENSP",
                     "SWISSPROT", "TREMBL", "UNIPARC", "SOURCE", "GENE_PHENO",
                     "SIFT", "PolyPhen", "DOMAINS", "miRNA", "HGVS_OFFSET",
                     "AF", "AFR_AF", "AMR_AF", "EAS_AF", "EUR_AF", "SAS_AF",
                     "AA_AF", "EA_AF", "gnomAD_AF", "gnomAD_AFR_AF",
                     "gnomAD_AMR_AF", "gnomAD_ASJ_AF", "gnomAD_EAS_AF",
                     "gnomAD_FIN_AF", "gnomAD_NFE_AF", "gnomAD_OTH_AF",
                     "gnomAD_SAS_AF", "MAX_AF", "MAX_AF_POPS", "CLIN_SIG",
                     "SOMATIC", "PHENO", "PUBMED", "MOTIF_NAME", "MOTIF_POS",
                     "HIGH_INF_POS", "MOTIF_SCORE_CHANGE", "CADD_PHRED",
                     "CADD_RAW", "REVEL", "LoFtool", "gnomADg", "gnomADg_AF",
                     "gnomADg_AF_afr", "gnomADg_AF_amr", "gnomADg_AF_eas",
                     "gnomADg_AF_fin", "gnomADg_AF_nfe", "gnomADg_AF_oth",
                     "gnomADg_popmax", "gnomADg_AF_popmax", "GERP_RS",
                     "phastCons100way_vertebrate", "phyloP100way_vertebrate",
                     "ClinVar", "ClinVar_CLNSIG", "ClinVar_CLNREVSTAT",
                     "ClinVar_CLNDN", "GroupID"])

    return variant_calls_dataframe


### Function 6:
    
def import_gene_list_annotations(path_to_annotation_data_folder, 
                                   list_of_annotation_file_names):
    """
    Reads a tab-delimited file with extra annotation information. The content
    of the files will always be two columns, where the first column is a list
    of gene names, and the second column the annotation information. 

    Parameters
    ----------
    path_to_annotation_data_folder : str
        The path to the annotation data folder
    
    list_of_annotation_file_names : list
        A list of strings, containing the anotation file names

    Returns
    -------
    list_of_annotation_file_dataframes : list
        A list of pandas dataframe objects
    """
                
    list_of_annotation_file_dataframes = []
        
    for annotation_file in list_of_annotation_file_names:
        path_to_ann_file = str(path_to_annotation_data_folder + 
                               str(annotation_file) + '.txt')
        ann_file_dataframe = pd.read_csv(path_to_ann_file, 
                                         sep = '\t', 
                                         header = 0)
        list_of_annotation_file_dataframes.append(ann_file_dataframe)
    
    return list_of_annotation_file_dataframes

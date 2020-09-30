"""
### File name:                                                              
process_split_vep_module.py                              

### General description of the script:                                      
This module contains functions for the following: Transformation of specific 
data columns, creation of new data columns based on variables contained in the 
original columns and based on data from the extra files.                                                                    
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Imports                                                                     #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# The script uses pandas for storing data in a dataframe object, and to perform
# various basic data analysis steps (extraction, subsetting etc.). 

import pandas as pd
import numpy as np

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Functions                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# Function 0:

def add_patient_ID_column(variant_calls_dataframe, ID):
    """
    Adds a new ID column to the variant_calls_dataframe.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    ID : str
        The sample ID to be added

    Returns
    -------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
    """
    
    variant_calls_dataframe['Patient_ID'] = ID
    
    return variant_calls_dataframe


# Function 1:
    
def transform_transcript_specific_columns(variant_calls_dataframe):
    """
    Transforms all transcript-specific columns from type string to list.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    variant_calls_transformed : pandas dataframe object
        The dataframe containing all variant calls
    """
    
    # Split the variant calls dataframe into two dataframes:
    # The first dataframe contains all the non-transcript specific columns
    # which are the first 11 columns. The second dataframe contains the tran-
    # script specific columns, which are the 89 remaining columns.
    non_transcript_specific_columns = variant_calls_dataframe.iloc[:, :11]
    transcript_specific_columns = variant_calls_dataframe.iloc[:, 11:99]
    groupID_column = variant_calls_dataframe.iloc[:, 99:]

    # Function to split strings into lists of strings
    def string_to_list(x):
        '''text'''
        return x.split(',')
    
    # Apply the function to all of the columns in the transcript specific data-
    # frame.
    transcript_specific_columns = transcript_specific_columns.applymap(
            string_to_list)

    # Merge the two dataframes again
    variant_calls_transformed = pd.concat([non_transcript_specific_columns, 
                                           transcript_specific_columns, 
                                           groupID_column], axis = 1, 
                                          sort = False)

    return variant_calls_transformed


# Function 2:

def annotate_variant_highest_impact_level(variant_calls_dataframe):
    """
    Identifies the highest impact level of each variant and adds a number 
    of extra columns with information on the primary impact level and 
    the indicies of transcripts with this impact level. Also adds new columns
    with the primary CADD Phred score and the primary gnomAD_AF.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    variant_calls_transformed : pandas dataframe object
        The dataframe containing all variant calls
    """
    
    # Function to identify the highest impact level of a given variant
    def extract_impact_level(impact_list):
        '''text'''
        if 'HIGH' in impact_list:
            impact = 'H'
            indices = [i for i, x in enumerate(impact_list) if x == "HIGH"]        
        elif 'MODERATE' in impact_list:
            impact = 'M'
            indices = [i for i, x in enumerate(impact_list) if x == "MODERATE"]
        elif 'LOW' in impact_list:
            impact = 'L'
            indices = [i for i, x in enumerate(impact_list) if x == "LOW"]
        elif 'MODIFIER' in impact_list:
            impact = 'O'
            indices = [i for i, x in enumerate(impact_list) if x == "MODIFIER"]
        return impact, indices
    
    # Apply the function 
    list_of_impacts = variant_calls_dataframe['Impact'].apply(
            extract_impact_level)

    # Convert to dataframe and rename columns
    dataframe_of_impacts = list_of_impacts.apply(pd.Series)
    dataframe_of_impacts.columns = ['Primary_impact', 'Impact_indicies']
    
    # Merge the two dataframes 
    variant_calls_transformed = pd.concat([variant_calls_dataframe, 
                                           dataframe_of_impacts], axis = 1, 
                                          sort = False)
    
    # Add column with CADD score
    list_of_CADD = []
    
    # Iterate over each row 
    for index, rows in variant_calls_transformed.iterrows(): 
        # Create list for the current row 
        cadd_row = [rows.CADD_PHRED[rows.Impact_indicies[0]]] 
        list_of_CADD.append(cadd_row)

    variant_calls_transformed['Primary_CADD_PHRED'] = list_of_CADD
    
    # Unlist CADD scores
    variant_calls_transformed['Primary_CADD_PHRED'] = \
        variant_calls_transformed['Primary_CADD_PHRED'].str[0]
    
    # Convert from string to float (. will be transformed to nan's)
    variant_calls_transformed['Primary_CADD_PHRED'] = \
        pd.to_numeric(variant_calls_transformed['Primary_CADD_PHRED'], 
                      errors = 'coerce')
    
    # Add column with gnom_AD allele frequencies
    list_of_gnomAD = []
    
    # Iterate over each row 
    for index, rows in variant_calls_transformed.iterrows(): 
        # Create list for the current row 
        gnomAD_row = [rows.gnomAD_AF[rows.Impact_indicies[0]]] 
        list_of_gnomAD.append(gnomAD_row)

    variant_calls_transformed['Primary_gnomAD_AF'] = list_of_gnomAD
    
    # Unlist CADD scores
    variant_calls_transformed['Primary_gnomAD_AF'] = \
        variant_calls_transformed['Primary_gnomAD_AF'].str[0]
    
    # Convert from string to float (. will be transformed to nan's)
    variant_calls_transformed['Primary_gnomAD_AF'] = \
        pd.to_numeric(variant_calls_transformed['Primary_gnomAD_AF'], 
                      errors = 'coerce')
    
    # Return the dataframe
    return variant_calls_transformed


# Function 3:

def annotate_variant_primary_consequence(variant_calls_dataframe):
    """
    Identifies the primary consequence for each variant and adds a number 
    of extra columns with information on the primary consequence  and 
    the indicies of transcripts with this consequence.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    variant_calls_transformed : pandas dataframe object
        The dataframe containing all variant calls
    """

    list_of_consequences = []
    list_of_indexes = []

    for i, row in variant_calls_dataframe.iterrows():
        # Get row consequences
        row_cons = row['Consequence']
        row_biotype = row['Biotype']
        # Check row impact level
        if row['Primary_impact'] == 'H':
            if any("transcript_ablation" in s for s in row_cons):
                list_of_consequences.append('transcript_ablation')
                index = [idx for idx, s in enumerate(row_cons) if 'transcript_ablation' in s][0]
                list_of_indexes.append(index)
            elif any("splice_acceptor_variant" in s for s in row_cons):
                list_of_consequences.append('splice_acceptor_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'splice_acceptor_variant' in s][0]
                list_of_indexes.append(index)
            elif any("splice_donor_variant" in s for s in row_cons):
                list_of_consequences.append('splice_donor_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'splice_donor_variant' in s][0]
                list_of_indexes.append(index)
            elif any("stop_gained" in s for s in row_cons):
                list_of_consequences.append('stop_gained')
                index = [idx for idx, s in enumerate(row_cons) if 'stop_gained' in s][0]
                list_of_indexes.append(index)
            elif any("frameshift_variant" in s for s in row_cons):
                list_of_consequences.append('frameshift_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'frameshift_variant' in s][0]
                list_of_indexes.append(index)
            elif any("stop_lost" in s for s in row_cons):
                list_of_consequences.append('stop_lost')
                index = [idx for idx, s in enumerate(row_cons) if 'stop_lost' in s][0]
                list_of_indexes.append(index)
            elif any("start_lost" in s for s in row_cons):
                list_of_consequences.append('start_lost')
                index = [idx for idx, s in enumerate(row_cons) if 'start_lost' in s][0]
                list_of_indexes.append(index)
            elif any("transcript_amplification" in s for s in row_cons):
                list_of_consequences.append('transcript_amplification')
                index = [idx for idx, s in enumerate(row_cons) if 'transcript_amplification' in s][0]
                list_of_indexes.append(index)
        elif row['Primary_impact'] == 'M':
            if any("inframe_insertion" in s for s in row_cons):
                list_of_consequences.append('inframe_insertion')
                index = [idx for idx, s in enumerate(row_cons) if 'inframe_insertion' in s][0]
                list_of_indexes.append(index)
            elif any("inframe_deletion" in s for s in row_cons):
                list_of_consequences.append('inframe_deletion')
                index = [idx for idx, s in enumerate(row_cons) if 'inframe_deletion' in s][0]
                list_of_indexes.append(index)
            elif any("missense_variant" in s for s in row_cons):
                list_of_consequences.append('missense_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'missense_variant' in s][0]
                list_of_indexes.append(index)
            elif any("protein_altering_variant" in s for s in row_cons):
                list_of_consequences.append('protein_altering_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'protein_altering_variant' in s][0]
                list_of_indexes.append(index)
            elif any("regulatory_region_ablation" in s for s in row_cons):
                list_of_consequences.append('regulatory_region_ablation')
                index = [idx for idx, s in enumerate(row_cons) if 'regulatory_region_ablation' in s][0]
                list_of_indexes.append(index)
        elif row['Primary_impact'] == 'L':
            if any("splice_region_variant" in s for s in row_cons):
                list_of_consequences.append('splice_region_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'splice_region_variant' in s][0]
                list_of_indexes.append(index)
            elif any("incomplete_terminal_codon_variant" in s for s in row_cons):
                list_of_consequences.append('incomplete_terminal_codon_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'incomplete_terminal_codon_variant' in s][0]
                list_of_indexes.append(index)
            elif any("start_retained_variant" in s for s in row_cons):
                list_of_consequences.append('start_retained_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'start_retained_variant' in s][0]
                list_of_indexes.append(index)
            elif any("stop_retained_variant" in s for s in row_cons):
                list_of_consequences.append('stop_retained_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'stop_retained_variant' in s][0]
                list_of_indexes.append(index)
            elif any("synonymous_variant" in s for s in row_cons):
                list_of_consequences.append('synonymous_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'synonymous_variant' in s][0]
                list_of_indexes.append(index)
        elif row['Primary_impact'] == 'O':
            if any("5_prime_UTR_variant" in s for s in row_cons):
                list_of_consequences.append('5_prime_UTR_variant')
                index = [idx for idx, s in enumerate(row_cons) if '5_prime_UTR_variant' in s][0]
                list_of_indexes.append(index)
            elif any("3_prime_UTR_variant" in s for s in row_cons):
                list_of_consequences.append('3_prime_UTR_variant')
                index = [idx for idx, s in enumerate(row_cons) if '3_prime_UTR_variant' in s][0]
                list_of_indexes.append(index)
            elif any("non_coding_transcript_exon_variant" in s for s in row_cons):
                list_of_consequences.append('non_coding_transcript_exon_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'non_coding_transcript_exon_variant' in s][0]
                list_of_indexes.append(index)
            elif any("non_coding_transcript_variant" in s for s in row_cons):
                list_of_consequences.append('non_coding_transcript_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'non_coding_transcript_variant' in s][0]
                list_of_indexes.append(index)
            elif any("regulatory_region_variant" in s for s in row_cons):
                if any("promoter" in s for s in row_biotype):
                    list_of_consequences.append('promoter_variant')
                    index = [idx for idx, s in enumerate(row_biotype) if 'promoter' in s][0]
                    list_of_indexes.append(index)
                elif any("TF_binding_site" in s for s in row_biotype):
                    list_of_consequences.append('TF_binding_site_variant')
                    index = [idx for idx, s in enumerate(row_biotype) if 'TF_binding_site' in s][0]
                    list_of_indexes.append(index)
                elif any("promotor_flanking_region" in s for s in row_biotype):
                    list_of_consequences.append('promotor_flanking_region_variant')
                    index = [idx for idx, s in enumerate(row_biotype) if 'promotor_flanking_region' in s][0]
                    list_of_indexes.append(index)
                elif any("CTCF_binding_site" in s for s in row_biotype):
                    list_of_consequences.append('CTCF_binding_site_variant')
                    index = [idx for idx, s in enumerate(row_biotype) if 'CTCF_binding_site' in s][0]
                    list_of_indexes.append(index)
                elif any("enhancer" in s for s in row_biotype):
                    list_of_consequences.append('enhancer_variant')
                    index = [idx for idx, s in enumerate(row_biotype) if 'enhancer' in s][0]
                    list_of_indexes.append(index)
                elif any("open_chromatin_region" in s for s in row_biotype):
                    list_of_consequences.append('open_chromatin_region')
                    index = [idx for idx, s in enumerate(row_biotype) if 'open_chromatin_region' in s][0]
                    list_of_indexes.append(index)
            elif any("TFBS_ablation" in s for s in row_cons):
                list_of_consequences.append('TFBS_ablation')
                index = [idx for idx, s in enumerate(row_cons) if 'TFBS_ablation' in s][0]
                list_of_indexes.append(index)
            elif any("TFBS_amplification" in s for s in row_cons):
                list_of_consequences.append('TFBS_amplification')
                index = [idx for idx, s in enumerate(row_cons) if 'TFBS_amplification' in s][0]
                list_of_indexes.append(index)
            elif any("TF_binding_site_variant" in s for s in row_cons):
                list_of_consequences.append('TF_binding_site_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'TF_binding_site_variant' in s][0]
                list_of_indexes.append(index)
            elif any("regulatory_region_amplification" in s for s in row_cons):
                list_of_consequences.append('regulatory_region_amplification')
                index = [idx for idx, s in enumerate(row_cons) if 'regulatory_region_amplification' in s][0]
                list_of_indexes.append(index)
            elif any("feature_elongation" in s for s in row_cons):
                list_of_consequences.append('feature_elongation')
                index = [idx for idx, s in enumerate(row_cons) if 'feature_elongation' in s][0]
                list_of_indexes.append(index)
            elif any("feature_truncation" in s for s in row_cons):
                list_of_consequences.append('feature_truncation')
                index = [idx for idx, s in enumerate(row_cons) if 'feature_truncation' in s][0]
                list_of_indexes.append(index)
            elif any("mature_miRNA_variant" in s for s in row_cons):
                list_of_consequences.append('mature_miRNA_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'mature_miRNA_variant' in s][0]
                list_of_indexes.append(index)
            elif any("NMD_transcript_variant" in s for s in row_cons):
                list_of_consequences.append('NMD_transcript_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'NMD_transcript_variant' in s][0]
                list_of_indexes.append(index)
            elif any("coding_sequence_variant" in s for s in row_cons):
                list_of_consequences.append('coding_sequence_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'coding_sequence_variant' in s][0]
                list_of_indexes.append(index)
            elif any("intergenic_variant" in s for s in row_cons):
                list_of_consequences.append('intergenic_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'intergenic_variant' in s][0]
                list_of_indexes.append(index)
            elif any("upstream_gene_variant" in s for s in row_cons):
                list_of_consequences.append('upstream_gene_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'upstream_gene_variant' in s][0]
                list_of_indexes.append(index)
            elif any("downstream_gene_variant" in s for s in row_cons):
                list_of_consequences.append('downstream_gene_variant')
                index = [idx for idx, s in enumerate(row_cons) if 'downstream_gene_variant' in s][0]
                list_of_indexes.append(index)

    # concatenate into new dataframe
    variant_calls_primary = pd.DataFrame(list(
            zip(list_of_consequences,list_of_indexes)), 
            columns =['primary_consequence', 'primary_index'])

    # Add to variant_call_dataframe
    variant_calls_transformed = pd.concat([variant_calls_dataframe, 
                                           variant_calls_primary], axis = 1, 
                                          sort = False)

    return variant_calls_transformed


# Function 4:

def add_transcript_number_column(variant_calls_dataframe):
    """
    Adds a new transcript number column to the variant_calls_dataframe.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
    """
        
    # Save as variable how many annotation groups each variant have
    # using the length of symbols list
    n_transcripts_per_var = \
     variant_calls_dataframe['Symbol'].apply(len)

    # Add new column with the total number of annotation group for each variant
    variant_calls_dataframe['number_of_transcripts'] = n_transcripts_per_var

    # Return the transformed variant call dataframe
    return variant_calls_dataframe


# Function 5:

def transform_groupID_column(variant_calls_dataframe):
    """
    Transforms the groupID column from type string to a list of strings. This
    list is further transformed into a fixed number of new columns. Each column 
    corresponds to a groupID.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    variant_calls_transformed : pandas dataframe object
        The dataframe containing all variant calls
    """
 
    # Convert all & to ,
    variant_calls_dataframe["GroupID"] = \
        variant_calls_dataframe.GroupID.apply(lambda x: x.replace('&', ','))

    # Split by ,
    variant_calls_dataframe["GroupID"] = \
        variant_calls_dataframe.GroupID.apply(lambda x: x.split(','))

    # Remove duplicates
    def remove_duplicates(list_of_groupIDs):
        '''text'''
        res = [] 
        [res.append(x) for x in list_of_groupIDs if x not in res]
        return res
    
    # apply function
    variant_calls_dataframe["GroupID"] = \
        variant_calls_dataframe.GroupID.apply(remove_duplicates)

    # Save as variable how many annotation groups each variant have
    n_ann_group_per_var = variant_calls_dataframe['GroupID'].apply(len)

    # Split GroupID column into multiple columns
    variant_calls_groups = pd.DataFrame(
        variant_calls_dataframe.GroupID.values.tolist()).add_prefix("Group_")

    # Extract the list of gene symbols for the first annotation group for 
    # each variant
    gene_symbol_list = variant_calls_groups["Group_0"]

    # Split by _ and keep only the first part of the string (the gene symbol)
    gene_symbol_list = gene_symbol_list.apply(lambda x: x.split("_")).str[0]

    ### Remove gene symbol from all columns
    
    # Get a list of the group columns
    group_columns = list(variant_calls_groups) 

    # Iterate over the group columns and remove the gene symbol
    for i in group_columns:
        variant_calls_groups[i] = variant_calls_groups[i].str.split('_').str[1]

    # If the length of the variant_calls_groups dataframe is less than 13, (the
    # max possible number of annotation groups) not all annotation groups 
    # appear. To make concatenation of dataframes from  different patients 
    # easier, lacking group columns are inserted
    
    while variant_calls_groups.shape[1] != 13:
        column_num = variant_calls_groups.shape[1]
        column_name = str('Group_' + str(column_num))
        variant_calls_groups[column_name] = ""
        variant_calls_groups[column_name] = np.nan

    # Replace Python object type 'None' with NaN instead
    variant_calls_groups = variant_calls_groups.fillna(value = np.nan)

    ### Replace columns and add columns

    # Remove old groupID column in main dataframe and replace with new columns
    variant_calls_dataframe = \
        variant_calls_dataframe.drop(["GroupID"], axis = 1 )

    variant_calls_transformed = pd.concat(
            [variant_calls_dataframe, variant_calls_groups], axis=1, sort=False)
    
    # Add new column with the annotated gene symbol
    variant_calls_transformed['Annotated_gene_symbol'] = gene_symbol_list
        
    # Add new column with the primary annnotation group
    variant_calls_transformed['Primary_annotation_group'] = \
        variant_calls_transformed['Group_0']
    
    variant_calls_transformed['Primary_annotation_group'] = \
        variant_calls_transformed['Primary_annotation_group'].str.split(
                '.').str[0]

    # Add new column with the total number of annotation group for each variant
    variant_calls_transformed['Total_N_ann_groups'] = n_ann_group_per_var

    return variant_calls_transformed


# Function 6:

def annotate_with_gene_list_info(variant_calls_dataframe, 
                                 list_of_annotation_info_dataframes):
    """
    Add new columns to genes with extra annotation information. 

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    list_of_annotation_file_dataframes : list
        A list of pandas dataframe objects

        
    Returns
    -------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
    """
    curated_gene_list_annotation = list_of_annotation_info_dataframes[0]
    reactome_gene_list_annotation = list_of_annotation_info_dataframes[1]
    go_gene_list_annotation = list_of_annotation_info_dataframes[2]
    kegg_gene_list_annotation = list_of_annotation_info_dataframes[3]
    hpo_gene_list_annotation = list_of_annotation_info_dataframes[4]
    gado_gene_list_annotation = list_of_annotation_info_dataframes[5] 
    findzebra_gene_list_annotation = list_of_annotation_info_dataframes[6]
    string_gene_list_interaction_partners = list_of_annotation_info_dataframes[7]
    NCBI_entrez_id = list_of_annotation_info_dataframes[8]
    Uniprot_identifier = list_of_annotation_info_dataframes[9]

    # Add gene list annotation to the genes
    variant_calls_dataframe['ThromboGenomics_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            curated_gene_list_annotation.set_index('gene_symbol')['thrombogenomics_annotation'])
    
    variant_calls_dataframe['ISTH_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            curated_gene_list_annotation.set_index('gene_symbol')['ISTH_annotation'])
    
    reactome_gene_list_annotation = reactome_gene_list_annotation.groupby('gene_symbol')['reactome_annotation'].apply('\n'.join).reset_index()

    variant_calls_dataframe['Reactome_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            reactome_gene_list_annotation.set_index('gene_symbol')['reactome_annotation'])

    go_gene_list_annotation = go_gene_list_annotation.groupby('gene_symbol')['go_annotation'].apply('\n'.join).reset_index()

    variant_calls_dataframe['GO_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            go_gene_list_annotation.set_index('gene_symbol')['go_annotation'])

    kegg_gene_list_annotation = kegg_gene_list_annotation.groupby('gene_symbol')['kegg_annotation'].apply('\n'.join).reset_index()

    variant_calls_dataframe['KEGG_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            kegg_gene_list_annotation.set_index('gene_symbol')['kegg_annotation'])

    hpo_gene_list_annotation = hpo_gene_list_annotation.groupby('gene_symbol')['hpo_annotation'].apply('\n'.join).reset_index()

    variant_calls_dataframe['HPO_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            hpo_gene_list_annotation.set_index('gene_symbol')['hpo_annotation'])
    
    gado_gene_list_annotation = gado_gene_list_annotation.groupby('gene_symbol')['gado_annotation'].apply('\n'.join).reset_index()

    variant_calls_dataframe['GADO_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            gado_gene_list_annotation.set_index('gene_symbol')['gado_annotation'])
    
    findzebra_gene_list_annotation = findzebra_gene_list_annotation.groupby('gene_symbol')['findzebra_annotation'].apply('\n'.join).reset_index()
    
    variant_calls_dataframe['FindZebra_annotation'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            findzebra_gene_list_annotation.set_index('gene_symbol')['findzebra_annotation'])

    string_gene_list_interaction_partners = string_gene_list_interaction_partners.groupby('gene_symbol')['interaction_partner'].apply(' & '.join).reset_index()
    
    variant_calls_dataframe['String_interaction_partner'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            string_gene_list_interaction_partners.set_index('gene_symbol')['interaction_partner'])

    NCBI_entrez_id['Entrez_ID'] = NCBI_entrez_id['Entrez_ID'].apply(str)

    NCBI_entrez_id = NCBI_entrez_id.groupby('Gene_symbol')['Entrez_ID'].apply(' & '.join).reset_index()
    
    variant_calls_dataframe['NCBI_entrez_id'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            NCBI_entrez_id.set_index('Gene_symbol')['Entrez_ID'])

    Uniprot_identifier = Uniprot_identifier.groupby('Gene_symbol')['Uniprot_identifier'].apply(' & '.join).reset_index()
    
    variant_calls_dataframe['Uniprot_identifier'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            Uniprot_identifier.set_index('Gene_symbol')['Uniprot_identifier'])

    # Split by ,
    string_gene_list_interaction_partners["interaction_partner_list"] = \
        string_gene_list_interaction_partners.interaction_partner.apply(lambda x: x.split('&'))

    # Save as variable how many string interactions each gene have
    n_string_group_per_gene = string_gene_list_interaction_partners['interaction_partner_list'].apply(len)
    
    # Add new column with the total string interactions for each gene 
    string_gene_list_interaction_partners['Total_N_string_interactions'] = n_string_group_per_gene

    variant_calls_dataframe['Total_N_string_interactions'] = variant_calls_dataframe['Annotated_gene_symbol'].map(
            string_gene_list_interaction_partners.set_index('gene_symbol')['Total_N_string_interactions'])
    
    values = {'ThromboGenomics_annotation': 'No entry', 
              'ISTH_annotation': 'No entry', 
              'Reactome_annotation': 'No entry', 
              'GO_annotation': 'No entry', 
              'KEGG_annotation': 'No entry', 
              'HPO_annotation': 'No entry', 
              'GADO_annotation': 'No entry', 
              'FindZebra_annotation': 'No entry', 
              'String_interaction_partner': 'No entry'}
    
    variant_calls_dataframe = variant_calls_dataframe.fillna(value=values)
    
    return variant_calls_dataframe

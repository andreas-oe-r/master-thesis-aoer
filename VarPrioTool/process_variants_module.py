"""
### File name:                                                              
process_variants_module.py                              

### General description of the script:                                      
This module has functions for the prioritisation of variants, and the creation
of data structures that are used as input for the functions defined in 
create_report_module.py                                                        
                                                                          
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Imports                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# Description:
# The script uses pandas for storing data in a dataframe object, and to perform
# various basic data analysis steps (extraction, subsetting etc.). Some 
# Reportlab data structures are also intialised here.

import pandas as pd
import numpy as np
from reportlab.platypus import (Paragraph)
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.enums import TA_JUSTIFY

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Functions                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# Function 0: ...
    
def extract_variant_classes(variant_calls_dataframe):
    """
    This function extracts specific variant groups.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    Returns
    -------
    homo_low, homo_moderate, homo_high, hetero_low
    hetero_moderate, hetero_high, pos_compound, 
    primary_reg_vars, other_vars : pandas dataframe objects
        Dataframes containing specific variant groups
    """
    # Get all homozygote variants
    homo_var = variant_calls_dataframe[
            (variant_calls_dataframe.Genotype == "1/1") | 
            (variant_calls_dataframe.Genotype == "1|1") | 
            (variant_calls_dataframe.Genotype == "2/2") | 
            (variant_calls_dataframe.Genotype == "2|2")]

    # Get all homozygotes variants with predicted low impact
    homo_low = homo_var[(homo_var.Primary_impact == "L")]
        
    # Get all homozygotes variants with predicted moderate impact
    homo_moderate = homo_var[(homo_var.Primary_impact == "M")]
        
    # Get all homozygotes variants with predicted high impact
    homo_high = homo_var[(homo_var.Primary_impact == "H")]
        
    # Get all heterozygote variants
    hetero_var = variant_calls_dataframe[
            (variant_calls_dataframe.Genotype == "0/1") | 
            (variant_calls_dataframe.Genotype == "1/0") |
            (variant_calls_dataframe.Genotype == "0|1") |
            (variant_calls_dataframe.Genotype == "1|0") |
            (variant_calls_dataframe.Genotype == "1/2") |
            (variant_calls_dataframe.Genotype == "2/1") |
            (variant_calls_dataframe.Genotype == "1|2") |
            (variant_calls_dataframe.Genotype == "2|1")]
        
    # Get all heterozygote variants with predicted low impact
    hetero_low = hetero_var[(hetero_var.Primary_impact == "L")]
        
    # Get all heterozygote variants with prediced moderate impact
    hetero_moderate = hetero_var[(hetero_var.Primary_impact == "M")]
        
    # Get all heterozygote variants with prediced high impact
    hetero_high = hetero_var[(hetero_var.Primary_impact == "H")]
    
    # Get all possible compound heterozygote variants 
    pos_compound = variant_calls_dataframe.groupby(
            'Annotated_gene_symbol').filter(
            lambda g: len(g) > 1).drop_duplicates(
            subset=['Annotated_gene_symbol', 'Chromosome_Position'], 
            keep="first")

    # Get all variants with impact = "MODIFIER"
    other_vars = variant_calls_dataframe[
            (variant_calls_dataframe.Primary_impact == "O")]
    
    # Get 5' UTR, and 3' UTR variants for the gene ANKRD26 and 
    # non_coding_transcript_(exon) variants for the gene RNU4ATAC
    primary_reg_vars = variant_calls_dataframe[
    (variant_calls_dataframe.primary_consequence == "5_prime_UTR_variant") | 
    (variant_calls_dataframe.primary_consequence == "3_prime_UTR_variant") |
    (variant_calls_dataframe.primary_consequence == \
     "non_coding_transcript_exon_variant")]
    
    primary_reg_vars = primary_reg_vars[
            (primary_reg_vars.Annotated_gene_symbol == "ANKRD26") |
            (primary_reg_vars.Annotated_gene_symbol == "RNU4ATAC")]
    
    return (homo_low, homo_moderate, homo_high, hetero_low, hetero_moderate,
            hetero_high, pos_compound, primary_reg_vars, other_vars)


# Function 1: ...
    
def extract_group_list_variants(variant_calls_dataframe, groupID):
    """
    This functions extracts the variants acording to groupID (grp1 is the 
    curated gene panel, while grp2-6 is the expanded gene list). Also performs
    quality control of the genotype, and further filtering of grp2-6 based on
    allele frequency. 

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    groupID : str
        A string that is either 'grp1' or 'grp2-6', indiciating which group
        should be extracted
        
    Returns
    -------
    grp_list_variants, grp_failed_qc, grp_passed_qc, 
    grp_SNPs_passed, grp_indels_passed, grp_moderate_impact, 
    grp_other_vars, grp_primary_reg_vars, grp_high_impact, 
    grp_moderate_impact_selected, grp_low_impact : pandas dataframe object
        Dataframes containing specific variant groups
    """
    
    if groupID == 'grp1':
        # Extract all variants with grp1 as primary group
        grp_list_variants = variant_calls_dataframe.loc[
                variant_calls_dataframe['Primary_annotation_group'] == 'grp1']
    
    if groupID == 'grp2-6':
        # Extract all variants with grp2 to grp6 as primary group
        grp_list_variants =  variant_calls_dataframe[
                variant_calls_dataframe['Primary_annotation_group'].isin(
                ['grp2', 'grp3', 'grp4', 'grp5', 'grp6'])]    
    
    # Extract all where genotype quality is above 50
    grp_passed_qc = grp_list_variants[
            grp_list_variants.Genotype_Quality_GQ > 49] 

    # Extract all where genotype quality is below 50
    grp_failed_qc = grp_list_variants[
            grp_list_variants.Genotype_Quality_GQ <= 49] 

    if groupID == 'grp2-6':
        # Drop all variants with allele frequency higher than 0.0001
        grp_passed_qc = grp_passed_qc[
            (grp_passed_qc.Primary_gnomAD_AF < 0.0001) |
            (grp_passed_qc.Primary_gnomAD_AF.isna())] 

    # Get the number of SNPS and indels passed
    grp_SNPs_passed = grp_passed_qc[grp_passed_qc['Type'].isin(
            ['SNP', 'SNP,OVERLAP'])]    

    grp_indels_passed = grp_passed_qc[grp_passed_qc['Type'].isin(
            ['INDEL', 'INDEL,OVERLAP'])]    

    # Get the various variant classes
    (grp_homo_low, grp_homo_moderate, grp_homo_high,
     grp_het_low, grp_het_moderate, grp_het_high, 
     grp_pos_comp, grp_primary_reg_vars, 
     grp_other_vars) = extract_variant_classes(grp_passed_qc)
        
    # Merge high impact candidates into one list
    grp_high_impact = pd.concat([grp_homo_high, grp_het_high])
    
    # Merge moderate impact candidates into one list
    grp_moderate_impact = pd.concat([grp_homo_moderate, 
                                     grp_het_moderate])

    # Merge low impact candidates into one list
    grp_low_impact = pd.concat([grp_homo_low, grp_het_low])

    ### Selection of moderate impact variants for group2-6
    # Drop all mitochondrial mutations
    # Select all inframe_insertion and inframe_deletions 
    # Select top 10 missense variants (by in-silico damage prediction)

    if groupID == 'grp2-6':
        # Drop all mitochondrial mutations
        grp_moderate_impact_no_MT = grp_moderate_impact[
                ~grp_moderate_impact.Chromosome_Position.str.contains("MT")]
   
        # Select any inframe insertion or deletions and save in seperate dataframe
        grp_moderate_impact_inframe = grp_moderate_impact_no_MT[
            (grp_moderate_impact_no_MT.primary_consequence == "inframe_insertion") |
            (grp_moderate_impact_no_MT.primary_consequence == "inframe_deletion")]
        
        # Select missense variants
        grp_moderate_impact_missense = grp_moderate_impact_no_MT[
            (grp_moderate_impact_no_MT.primary_consequence == "missense_variant")]
        
        # Sort acording to CADD Phred scale 
        grp_moderate_impact_missense = \
            grp_moderate_impact_missense.sort_values(
                    by = ['Primary_CADD_PHRED'], 
                    ascending = False)
        
        # Select top 10
        grp_moderate_impact_missense = grp_moderate_impact_missense.head(n = 10)
        
        # Concatenate the the two dataframes such that inframe variants
        # are listed before missense variants
        if grp_moderate_impact_inframe.empty:
            grp_moderate_impact_selected = grp_moderate_impact_missense
        else:
            grp_moderate_impact_selected = pd.concat(
                    [grp_moderate_impact_inframe, grp_moderate_impact_missense])

    if groupID == 'grp1':
        return (grp_list_variants, grp_failed_qc, grp_passed_qc, 
                grp_SNPs_passed, grp_indels_passed, grp_high_impact, 
                grp_moderate_impact, grp_low_impact, grp_primary_reg_vars,
                grp_other_vars)

    if groupID == 'grp2-6':
        return (grp_list_variants, grp_failed_qc, grp_passed_qc, 
                grp_SNPs_passed, grp_indels_passed, grp_moderate_impact, 
                grp_other_vars, grp_primary_reg_vars, grp_high_impact, 
                grp_moderate_impact_selected, grp_low_impact)


# Function 2: ...

def generate_summary_table(variant_calls_dataframe, groupID):
    """
    This function creates the summary table.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    groupID : str
        A string that is either 'grp1' or 'grp2-6'
        
    Returns
    -------
    summary_table : list object
        A list containing the fields to be inserted in a table in the 
        PDF report
    """
    if groupID == 'grp1':
#        summary_table_variants = pd.concat(variant_calls_dataframe[5:9])
        summary_table_variants = pd.concat([
                variant_calls_dataframe[5],
                variant_calls_dataframe[6],
                variant_calls_dataframe[7],
                variant_calls_dataframe[8]])
        
    if groupID == 'grp2-6':
#        summary_table_variants = pd.concat(variant_calls_dataframe[8:10])
        summary_table_variants = pd.concat([
                variant_calls_dataframe[8],
                variant_calls_dataframe[9]])
        
    # Add new column with numbering of the variants
    summary_table_variants['Variant_row_num'] = np.arange(
            len(summary_table_variants))

    summary_table = [['Num', 
                      'Chr:Pos',
                      'Gene', 
                      'Type', 
                      'Genotype',
                      'Worst consequence',
                      'GQ',
                      'DP',
                      'AD']]
  
    # Iterate over each row 
    for index, rows in summary_table_variants.iterrows(): 
        # Create list for the current row 
        variant_row = [str('Var_' + str(rows.Variant_row_num)),
                       rows.Chromosome_Position,
                       rows.Annotated_gene_symbol, 
                       rows.Type, 
                       rows.Genotype, 
                       rows.primary_consequence,
                       rows.Genotype_Quality_GQ,
                       rows.Depth_of_coverage_DP,
                       rows.Allele_depth_AD]
        summary_table.append(variant_row)

    return summary_table


# Function 3: ... 

def generate_stats_table(variant_calls_dataframe, groupID):
    """
    This function creates the statistics table.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    groupID : str
        A string that is either 'grp1' or 'grp2-6'
        
    Returns
    -------
    stats_table : list object
        A list containing the fields to be inserted in a table in the 
        PDF report
    """
    # Number of unique genes present
    n_genes = len(variant_calls_dataframe[0]['Annotated_gene_symbol'].unique())

    if groupID == 'grp1':
        genes_per = (n_genes / 126) * 100
    
    if groupID == 'grp2-6':
        genes_per = (n_genes / 2173) * 100

    # Variants before QC
    n_var = len(variant_calls_dataframe[0])
    
    # Variants after QC
    n_var_qc = len(variant_calls_dataframe[2])
    
    # Percentage remaning after QC
    var_per_qc = (len(variant_calls_dataframe[2]) / 
                 len(variant_calls_dataframe[0]) * 100)
    
    # Number of SNPs
    n_snps = len(variant_calls_dataframe[3])
    
    # Number of indels
    n_indels = len(variant_calls_dataframe[4])
    
    # Number of transcripts
    n_transcripts = variant_calls_dataframe[2]['number_of_transcripts'].sum()
    
    # Number of transcripts on average per variant
    transcripts_avg = n_transcripts / n_var_qc
    
    if groupID == 'grp1':
        group_genes_str = 'Curated gene list total:'
        group_genes = 125
        group_genes_pres = 'Curated gene list present:'
    
    if groupID == 'grp2-6':
        group_genes_str = 'Expanded gene list total:'
        group_genes = 2173
        group_genes_pres = 'Expanded gene list present:'
  
    stats_table = [[group_genes_str, group_genes],
                   [group_genes_pres, 
                str(str(n_genes) + ' (' + str(round(genes_per, 2)) + '%)')],
                   ['Variants total:', n_var],
                   ['Variants passed QC:', 
                str(str(n_var_qc) + ' (' + str(round(var_per_qc, 2)) + '%)')],
                    ['SNPs passed QC', n_snps],
                    ['Indels passed QC', n_indels],
                    ['Transcripts in total', n_transcripts],
                    ['Transcrips per variant (avg)', 
                     round(transcripts_avg, 2)]]
    
    return stats_table


# Function 4: ...

def generate_expanded_info_table(variant_calls_dataframe, groupID):
    """
    This function creates the expanded info table for each variant, and stores
    all tables in an overall list.

    Parameters
    ----------
    variant_calls_dataframe : pandas dataframe object
        The dataframe containing all variant calls
        
    groupID : str
        A string that is either 'grp1' or 'grp2-6'
        
    Returns
    -------
    stats_table : list object
        A list containing lists of the fields to be inserted in a table in the 
        PDF report
    """
    if groupID == 'grp1':
#        variant_table_variants = pd.concat(variant_calls_dataframe[5:9])
        variant_table_variants = pd.concat([
                variant_calls_dataframe[5],
                variant_calls_dataframe[6],
                variant_calls_dataframe[7],
                variant_calls_dataframe[8]])
                       
    if groupID == 'grp2-6':
#        variant_table_variants = pd.concat(variant_calls_dataframe[8:10])
        variant_table_variants = pd.concat([
                variant_calls_dataframe[8],
                variant_calls_dataframe[9]])
                       
    # Add new column with numbering of the variants
    variant_table_variants['Variant_row_num'] = np.arange(
            len(variant_table_variants))

    # Add new list to store the tables
    variant_tables_list = []

    # Define column header rows (row 1, row 3, row 5, row 7)
    
    table_1 = ['General variant information:']
    
    table_2_row_0 = ['Chr:Pos', 
                     'Gene symbol', 
                     'Type', 
                     'Genotype',
                     'Worst consequence', 
                     'GQ', 
                     'DP', 
                     'AD']
   
    table_x1_row_0 = ['Exon', 
                      'Intron',
                      'Variant class',
                      'Exisiting variation']

    table_x2_row_0  = ['Feature type', 
                       'Biotype',
                       'gnomAD AF', 
                       'Number of affected transcripts'] 
    
    table_x3_row_0  = ['CADD Phred',
                       'PolyPhen ',
                       'Sift',
                       'ClinVar significance',]
    
    table_x4_row_0  = ['PhastCons', 
                       'phyloP', 
                       'Ensembl gene', 
                       'Ensembl Feature']

    table_5 = ['Sequence information:']
        
    table_7 = ['Gene list specific annotation information:']
    
    table_9 = ['Gene list annotation links']
        
    # Set normal style
    styles = getSampleStyleSheet()
    styles.add(ParagraphStyle(name='normal', alignment=TA_JUSTIFY))
    styles["normal"].fontSize = 10 
    styles["normal"].leading = 10

    # Iterate over each row in the dataframe (corresponding to a variant)
    for index, rows in variant_table_variants.iterrows(): 
        # Create list for the current row  / variant
        variant_table = []
        # Create table_0 (empty, will be filled by header created during
        # generation of the report)
        table_0 = []
        # Append to 
        variant_table.append(table_0)
        variant_table.append(table_1)
        # Create table 2
        table_2 = []
        table_2.append(table_2_row_0)
        table_2_row_1 = [rows.Chromosome_Position,
                         rows.Annotated_gene_symbol, 
                         rows.Type, 
                         rows.Genotype, 
                         rows.primary_consequence,
                         rows.Genotype_Quality_GQ,
                         rows.Depth_of_coverage_DP,
                         rows.Allele_depth_AD]
        table_2.append(table_2_row_1)
        variant_table.append(table_2)
        # Create table x1
        table_x1 = []
        table_x1.append(table_x1_row_0)
        table_x1_row_1 = [rows.Exon[rows.Impact_indicies[0]],
                         rows.Intron[rows.Impact_indicies[0]],
                         rows.VARIANT_CLASS[rows.Impact_indicies[0]],
                         rows.Existing_variation[rows.Impact_indicies[0]]]
        table_x1.append(table_x1_row_1)
        variant_table.append(table_x1)
        # Create table x2
        table_x2 = []
        table_x2.append(table_x2_row_0)
        table_x2_row_1 = [rows.Feature_type[rows.Impact_indicies[0]],
                          rows.Biotype[rows.Impact_indicies[0]],
                          rows.gnomAD_AF[rows.Impact_indicies[0]],
                          rows.number_of_transcripts]
        table_x2.append(table_x2_row_1)
        variant_table.append(table_x2)
        # Create table x3
        table_x3 = []
        table_x3.append(table_x3_row_0)
        table_x3_row_1 = [str(rows.Primary_CADD_PHRED),
                          rows.PolyPhen[rows.Impact_indicies[0]],
                          rows.SIFT[rows.Impact_indicies[0]],
                          rows.ClinVar_CLNSIG[rows.Impact_indicies[0]]]
        table_x3.append(table_x3_row_1)
        variant_table.append(table_x3)
        # Create table x4
        table_x4 = []
        table_x4.append(table_x4_row_0)
        table_x4_row_1 = [rows.phastCons100way_vertebrate[rows.Impact_indicies[0]],
                          rows.phyloP100way_vertebrate[rows.Impact_indicies[0]],
                          rows.Gene[rows.Impact_indicies[0]],
                          rows.Feature[rows.Impact_indicies[0]]]
        table_x4.append(table_x4_row_1)
        variant_table.append(table_x4)
        # Create table 5
        variant_table.append(table_5)
        # Create table 6
        table_6 = [['Reference', rows.Reference],
                   ['Alternative', rows.Alternative],
                   ['Allele', rows.Allele[rows.Impact_indicies[0]]],
                   ['HGVSc', rows.HGVSc[rows.Impact_indicies[0]]],
                   ['HGVSp', rows.HGVSp[rows.Impact_indicies[0]]]]
        variant_table.append(table_6)
        # Create table 7
        variant_table.append(table_7)
        # Create table 8
        table_8 = []
        if rows.ThromboGenomics_annotation != 'No entry':
            table_8_row_2 = ['ThromboGenomics annotation:', 
                             rows.ThromboGenomics_annotation]
            table_8.append(table_8_row_2)
        if rows.ISTH_annotation != 'No entry':
            table_8_row_3 = ['ISTH annotation:', 
                             rows.ISTH_annotation]
            table_8.append(table_8_row_3)
        if rows.Reactome_annotation != 'No entry':
            table_8_row_4 = ['Reactome pathway annotation:', 
                             rows.Reactome_annotation]
            table_8.append(table_8_row_4)
        if rows.GO_annotation != 'No entry':
            table_8_row_5 = ['Gene Ontology annotation:', 
                             rows.GO_annotation]
            table_8.append(table_8_row_5)
        if rows.KEGG_annotation != 'No entry':
            table_8_row_6 = ['KEGG pathway annotation:', 
                             rows.KEGG_annotation]
            table_8.append(table_8_row_6)
        if rows.HPO_annotation != 'No entry':
            table_8_row_7 = ['Human Phenotype Ontology annotation:', 
                             rows.HPO_annotation]
            table_8.append(table_8_row_7)
        if rows.GADO_annotation != 'No entry':
            table_8_row_8 = ['GADO annotation:', 
                             rows.GADO_annotation]
            table_8.append(table_8_row_8)
        if rows.FindZebra_annotation != 'No entry':
            table_8_row_9 = ['FindZebra annotation:', 
                             rows.FindZebra_annotation]
            table_8.append(table_8_row_9)
        if rows.String_interaction_partner != 'No entry':
            table_8_row_10 = ['STRING interactions w. curated genes:', 
                              rows.String_interaction_partner]
            table_8.append(table_8_row_10)
        table_8_row_1 = ['Annotation sources total:', len(table_8)]
        table_8.append(table_8_row_1)
        variant_table.append(table_8)
        # Create table 9
        variant_table.append(table_9)
        # Create table 10
        table_10 = []
        # Create table 10 – row 1
        table_10_row_1 = ['NCBI gene information:', 
                          Paragraph(
                                  str(
                            '<a href="https://www.ncbi.nlm.nih.gov/gene?term='
                                      + rows.NCBI_entrez_id 
                                      + '''" color="blue">Link to gene entry at
                                      NCBI gene</a>'''),
                                      style = styles["Normal"])]
        table_10.append(table_10_row_1)
        # Create table 10 – row 2
        table_10_row_2 = ['Ensembl gene information:', 
                          Paragraph(
                                  str(
    '<a href="http://grch37.ensembl.org/Homo_sapiens/Gene/Summary?db=core;g='
                                      + rows.Gene[rows.Impact_indicies[0]]
                                      + '''" color="blue">Link to gene entry at
                                      Ensembl genes</a>'''),
                                      style = styles["Normal"])]
        table_10.append(table_10_row_2)
        # Create table 10 – row 3
        table_10_row_3 = ['Unitprot protein information:', 
                          Paragraph(
                                  str(
                            '<a href="https://www.uniprot.org/uniprot/'
                                      + rows.Uniprot_identifier
                                      + '''" color="blue">Link to protein entry 
                                      at Uniprot</a>'''),
                                      style = styles["Normal"])]
        table_10.append(table_10_row_3)
        # Append row tables to the overal table
        variant_table.append(table_10)
        # Create table 11
        table_11 = []
        # Create table 10 – row 1
        if rows.Reactome_annotation == 'No entry':
            hemostasis_ann = '(gene is not part of hemostasis pathway)'
        else:
            hemostasis_ann = '(gene is part of hemostasis pathway)'
        table_11_row_1 = ['Reactome gene/reactions entries:', 
                          Paragraph(
                                  str(
                            '<a href="https://reactome.org/content/query?q='
                                      + rows.Uniprot_identifier
                                      + '''&species=Homo+sapiens" 
                                      color="blue">Link to search result at
                                      Reactome ''' + hemostasis_ann + '</a>'),
                                      style = styles["Normal"])]
        table_11.append(table_11_row_1)
        # Create table 10 – row 2
        if rows.GO_annotation != 'No entry':
            table_11_row_2 = ['Gene Ontology hemostasis overview:', 
                              Paragraph(
                                      str(
                    '''<a href="https://www.ebi.ac.uk/QuickGO/term/GO:0007599" 
                    color="blue">Link to hemostasis gene ontology</a>'''),
                                      style = styles["Normal"])]
            table_11.append(table_11_row_2)
        # Create table 10 – row 3
        if rows.KEGG_annotation != 'No entry':
            web_adress_1 = 'https://www.genome.jp/kegg-bin/show_pathway?map=hsa04611&show_description=show'
            web_adress_2 = 'https://www.genome.jp/kegg-bin/show_pathway?hsa04610'
            if 'platelet_activation' and 'complement' in rows.KEGG_annotation:
                table_11_row_3 = ['KEGG pathway overview:', 
                                  Paragraph(
                                          str('<a href="' + web_adress_1 + 
                                              '''" color="blue">Link to KEGG 
                                              platelet pathway</a>''' + ' | ' + 
                                              '<a href="' + web_adress_2 + 
                                              '''" color="blue">Link to KEGG 
                                              complement pathway</a>''' ),
                                              style = styles["Normal"])]
                table_11.append(table_11_row_3)
            elif 'platelet_activation' in rows.KEGG_annotation:
                table_11_row_3 = ['KEGG pathway overview:', 
                                  Paragraph(
                                          str('<a href="' + web_adress_1 + 
                                              '''" color="blue">Link to KEGG 
                                              platelet pathway</a>'''),
                                              style = styles["Normal"])]
                table_11.append(table_11_row_3)
            elif 'complement' in rows.KEGG_annotation:
                table_11_row_3 = ['KEGG pathway overview:', 
                                  Paragraph(
                                          str('<a href="' + web_adress_2 + 
                                              '''" color="blue">Link to KEGG 
                                              complement pathway</a>'''),
                                              style = styles["Normal"])]
                table_11.append(table_11_row_3)
        # Create table 11 – row 4
        if rows.HPO_annotation != 'No entry':
            table_11_row_4 = ['Human Phenotype Ontology:', 
                              Paragraph(
                            str('<a href="https://hpo.jax.org/app/browse/gene/'
                                      + rows.NCBI_entrez_id 
                                      + '''" color="blue">Link to gene entry at
                                      Human Phenotype Ontology 
                                      </a>'''), 
                                      style = styles["Normal"])]
            table_11.append(table_11_row_4)
        # Create table 1 – row 5
        if rows.GADO_annotation != 'No entry':
            table_11_row_5 = ['GADO gene information:', 
                              Paragraph(
                            str('<a href="https://www.genenetwork.nl/gene/'
                                      + rows.Annotated_gene_symbol 
                                      + '''" color="blue">Link to gene entry at
                                      GADO</a>'''), 
                                      style = styles["Normal"])]
            table_11.append(table_11_row_5)
        # Create table 11 – row 6
        if rows.FindZebra_annotation != 'No entry':
            table_11_row_6 = ['FindZebra disease search results:', 
                              Paragraph(
                            str('<a href="https://www.findzebra.com/genes?q='
                                      + rows.Annotated_gene_symbol 
                                      + '''" color="blue">Link to gene entry at
                                      FindZebra</a>'''), 
                                      style = styles["Normal"])]
            table_11.append(table_11_row_6)
        # Create table 1 – row 7
        if rows.String_interaction_partner != 'No entry':
            table_11_row_7 = ['STRING protein-protein interactions:', 
                              Paragraph(
                            str('<a href="https://string-db.org/network/'
                                      + rows.ENSP[rows.Impact_indicies[0]]
                                      + '''" color="blue">Link to protein PPI
                                      network at STRING database</a>'''), 
                                      style = styles["Normal"])]
            table_11.append(table_11_row_7)
        # Append row tables to the overal table
        variant_table.append(table_11)
        # Append overall table to the general list
        variant_tables_list.append(variant_table)
        
    return variant_tables_list


# Function 8: ...

def generate_IGV_batch_file(variant_calls_curated, variant_calls_expanded,
                            path_to_IGV_batch_file, path_to_snapshot_directory,
                            path_to_cram_file):
    """
    This function creates the igv batch file. 

    Parameters
    ----------
        
    Returns
    -------
    """
    # Concatenate the prioritised variants of the curated gene panel
    IGV_var_curated = pd.concat(variant_calls_curated[5:9])

    # Concatenate the prioritised variants of the expanded gene list
    IGV_var_exp = pd.concat(variant_calls_expanded[8:10])
    
    # Concatenate the curated variants and the expanded variants into a single
    # dataframe
    list_of_IGV_var = [IGV_var_curated, IGV_var_exp]
    IGV_var = pd.concat(list_of_IGV_var)
    
    # Drop all but the first column
    IGV_var = IGV_var[['Chromosome_Position']]

    # Split into two columns
    IGV_var[['Chromosome','Position']] = \
    IGV_var.Chromosome_Position.str.split(':', expand=True)

    # Convert position column from type string to int
    IGV_var["Position"] = IGV_var["Position"].astype(str).astype(int)

    # Add one column with position minus 100 bp
    IGV_var['Position_minus_100_bp'] = IGV_var['Position'] - 100
    
    # Add one column with position plus 100 bp
    IGV_var['Position_plus_100_bp'] = IGV_var['Position'] + 100

    # Create IGV batch file
    IGV_batch_file = open(path_to_IGV_batch_file, 'w')
    
    # Set the inital lines of the file
    initial_lines = [
     'new\n', 
     'snapshotDirectory ' + path_to_snapshot_directory + '\n', 
     'load ' + path_to_cram_file + '\n', 
     'Genome hg19\n',
     'load path/to/igv/hg19.genome\n',
     'load path/to/igv/genomes/hg19/cytoBand.txt\n', 
     'load /path/to/igv/genomes/hg19/hg19_alias.tab\n',
     'load /path/to/igv/genomes/hg19/property.txt\n',
     'load path/to/igv/genomes/hg19/refGene.txt\n']
        
    # Write initial lines to file
    IGV_batch_file.writelines(initial_lines)
        
    # Set count variable to used as header (e.g. Var 0, Var 1, Var 2 ...)
    i = 0

    # Write lines for each variant to file
    for index, row in IGV_var.iterrows():
        # Write header
        var_id = i
        # Get the chromosome, and position-minus and postion-plus of the 
        # variant
        chrom = str('chr' + row['Chromosome'])
        pos_minus = str(row['Position_minus_100_bp'])
        pos_plus = str(row['Position_plus_100_bp'])
        # Write the position line to the file
        variant_lines = ['goto ' + chrom + ':' + pos_minus + '-' + pos_plus + 
                         '\n']
        IGV_batch_file.writelines(variant_lines)
        # Write the endlines for the variants
        end_lines = ['sort quality\n', 'squish\n', 'snapshot ' + str(var_id) + 
                     '.png\n']
        IGV_batch_file.writelines(end_lines)
        # Increment the header count variable 
        i = i + 1
    
    # Write exit lines to file
    IGV_batch_file.write('exit')
    
    # Close the file
    IGV_batch_file.close()
        
    return 0

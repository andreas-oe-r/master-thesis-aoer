"""
### File name:                                                              
generate_variant_prioritisation_report_main.py                              

### General description of the script:                                      
This python script contains the main function that performs a series of steps 
to generate a PDF-report with results, by analysing data from a variant call 
format (VCF) file. The data that is analysed is imported from a tab-delimited 
data-file  which has been produced using the BCFtools plugin split-vep, which 
takes as input a VCF file that has been annotated with Ensembl Variant Effect 
Predictor (VEP).

The script utilises a number of functions defined in four auxillary modules for
import, transformation, annotation and data-mining of the data-file (and 
auxillary files) and generation of a PDF report detaling the results. Important 
packages used in the auxillary modules include pandas, numpy and reportlab. In 
this script, the overall stepwise procedure is documented. For documentation of
the various functions, see the specific modules.                                                                                                                                          
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Import of modules                                                           #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

import import_files_module as impfil
import process_split_vep_module as provep
import process_variants_module as provar
import create_report_module as genrep

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Defining the main function                                                  #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

def main():
    """
    This is the entry point of the Python program, and the execution of the
    program consists of calling this function (see line xxx).
    """
        
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 0: Assigning sample-specific variables                             #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    ### Part 0.0: Get patient ID from file
    
    patient_ID = impfil.read_patient_ID_file(
            path_to_patient_ID_file = str('file.txt'))
    
    ### Part 0.1: Get the list of patient-specific project sample ID's from 
    # file
    project_sample_ID_list = impfil.read_project_sample_file(
            path_to_project_sample_file = str('file.txt'))

    ### Part 0.2: Get the appropriate project sample ID from the patient ID
    project_sample_ID = impfil.get_sample_ID_from_patient_ID(
            patient_ID = patient_ID,
            project_sample_ID_list = project_sample_ID_list)
        
    ### Part 0.3: Path to patient-specific folder and files
    # Assign path to patient-specific folder, variant analysis report, IGV 
    # screenshots file, the meta-analysis directory, and paths for the two 
    # export tables for the meta-analysis.
        
    path_to_report_folder = str('path/on/server/')

    # Assign path to the IGV snapshots directory
    path_to_IGV_snapshots_directory = str('path/on/server/')
        
    # Assign path to variant analysis report 
    path_to_variant_analysis_report = str(path_to_report_folder + 
                                          'variant_analysis_report_' +
                                           str(patient_ID) +
                                          '.pdf')
        
    # Assign path to the meta-analysis directory
    path_to_meta_analysis_directory = str('path/on/server/')

    curated_variants_export_table = str(path_to_meta_analysis_directory + 
                                        str(patient_ID) + '_curated.txt')
    
    expanded_variants_export_table = str(path_to_meta_analysis_directory + 
                                         str(patient_ID) + '_expanded.txt')
    
    ### Part 0.4: Assign path to VCF data
    # Path to VCF data
    path_to_split_vep_file = str('path/on/server/')
            
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 1: Importing data-file and annotation data                         #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    ### Part 1.0: Import of tab-delimited file with variant calls information

    variant_calls = impfil.import_split_vep_file(
            path_to_split_vep_file_str = path_to_split_vep_file)
        
    ### Part 1.1: Import annotation data from gene list sources
    
    # Create list with file names
    annotation_file_names_list = ['curated_gene_list_gene_annotations',
                                  'reactome_gene_annotations',
                                  'go_gene_annotations',
                                  'kegg_gene_annotations',
                                  'hpo_gene_annotations',
                                  'gado_gene_annotations',
                                  'findzebra_gene_annotations',
                                  'string_gene_interaction_partners',
                                  'NCBI_entrez_id',
                                  'uniprot_identifier_list']
    
    # Import all of the files and store results as a list of dataframes
    annotation_info_dataframes_list = impfil.import_gene_list_annotations(
            path_to_annotation_data_folder = str('path/on/server/'),
            list_of_annotation_file_names = annotation_file_names_list)
            
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 2: Transformation and processing of the variant call dataframe     #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    ### Part 2.0: Add new column with the patient ID for use in later 
    # meta-analysis of the patient cohorte
    variant_calls = provep.add_patient_ID_column(
            variant_calls_dataframe = variant_calls,
            ID = project_sample_ID)

    ### Part 2.1: Transformation of all transcript-specific columns (expect the
    # groupID column) from type string/float to list of strings/floats (such 
    # that each transcript-specific column has one entry for each transcript 
    # in a list)
    variant_calls = provep.transform_transcript_specific_columns(
            variant_calls_dataframe = variant_calls)

    ### Part 2.2: Add new columns for each variants with information on the
    # highest transcript impact level of the transcipts and the indicies of the 
    # transcripts with this level (one impact level for each variant and 
    # x-number of indicies, depending on the number of transcripts)
    variant_calls = provep.annotate_variant_highest_impact_level(
            variant_calls_dataframe = variant_calls)

    ### Part 2.3: Add new column with the primary (worst) consequence of each
    # variant, and the index of the transcript with this consequence
    variant_calls = provep.annotate_variant_primary_consequence(
            variant_calls_dataframe = variant_calls)

    ### Part 2.4:
    # Add new column with the number of transcripts present for each variant
    variant_calls = provep.add_transcript_number_column(
            variant_calls_dataframe = variant_calls)

    ### Part 2.5:
    # Transformation of the groupID column from one string of groupIDs to 
    # multiple columns, where each column represents one annotation group, and
    # addition of new column with the primary annotation group and new column
    # with the total number of annotation groups
    variant_calls = provep.transform_groupID_column(
            variant_calls_dataframe = variant_calls)
    
    ### Part 2.6:
    # Add gene list annotation to the genes
    variant_calls = provep.annotate_with_gene_list_info(
        variant_calls_dataframe = variant_calls,
        list_of_annotation_info_dataframes = annotation_info_dataframes_list)
        
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 3: Extraction and processing of the variants by group              #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    # Step 3.0: Processing of curated gene list variants (groupID = group 1)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Extraction of variants from curated gene list / gene panel (group 1)
    # and quality filtering based on genotype quality (all variants with 
    # genotype quality < 50 are not considered)
    curated_variants = provar.extract_group_list_variants(
                variant_calls_dataframe = variant_calls,
                groupID = 'grp1')
        
    # Generation of tables to be inserted into the report
    cur_summary_table = provar.generate_summary_table(
            variant_calls_dataframe = curated_variants,
            groupID = 'grp1')    
    
    cur_stats_table = provar.generate_stats_table(
            variant_calls_dataframe = curated_variants,
            groupID = 'grp1')    
    
    cur_var_table = provar.generate_expanded_info_table(
            variant_calls_dataframe = curated_variants,
            groupID = 'grp1')       
            
    # Step 3.1: Processing of expanded gene list variants (groupID = group 2-6)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Extraction of variants from group 2-6 (expanded gene list)
    # and quality filtering based on genotype quality (all variants with 
    # genotype quality < 50 are not considered)
    expanded_variants = provar.extract_group_list_variants(
                variant_calls_dataframe = variant_calls,
                groupID = 'grp2-6')
            
    # Generation of tables to be inserted into the report
    exp_summary_table = provar.generate_summary_table(
            variant_calls_dataframe = expanded_variants,
            groupID = 'grp2-6')    
    
    exp_stats_table = provar.generate_stats_table(
            variant_calls_dataframe = expanded_variants,
            groupID = 'grp2-6')    
    
    exp_var_table = provar.generate_expanded_info_table(
            variant_calls_dataframe = expanded_variants,
            groupID = 'grp2-6')    

    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 4: Generate variant analysis report PDF                            #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    # Generate the PDF object
    variant_analysis_report = genrep.generate_pdf_object(
            path_to_pdf_file_str = path_to_variant_analysis_report)

    # Create a list of elements to included in the report
    list_of_report_elements = [cur_summary_table, cur_stats_table,
                               exp_summary_table, exp_stats_table, 
                               cur_var_table, exp_var_table]

    # Generate the report
    genrep.generate_variant_analysis_report(
            patient_id = project_sample_ID,
            pdf_object = variant_analysis_report,
            pdf_elements = list_of_report_elements,
            path_to_IGV_snapshots_directory = path_to_IGV_snapshots_directory)
    
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 5: Export variant table for meta analysis                          #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    curated_variants_dataframe = curated_variants[2]
    curated_variants_dataframe.to_csv(curated_variants_export_table, 
                                      sep = '\t',
                                      index = False)
    
    expanded_variants_dataframe = expanded_variants[2]
    expanded_variants_dataframe.to_csv(expanded_variants_export_table,
                                      sep = '\t',
                                      index = False)

    return 0

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Execution of the main function                                              #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

if __name__=="__main__":
    main()

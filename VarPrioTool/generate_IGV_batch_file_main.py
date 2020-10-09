"""
### File name:                                                              
generate_variant_prioritisation_report_main.py                              

### General description of the script:                                      
This python script contains the main function that performs a series of steps 
to generate an IGV batch file, by analysing data from a variant call 
format (VCF) file. The data that is analysed is imported from a tab-delimited 
data-file  which has been produced using the BCFtools plugin split-vep, which 
takes as input a VCF file that has been annotated with Ensembl Variant Effect 
Predictor (VEP).

The script utilises a number of functions defined in three auxillary modules 
for import, transformation, annotation and data-mining of the data-file (and 
auxillary files). Important packages used in the auxillary modules include 
pandas, numpy and reportlab. In this script, the overall stepwise procedure is 
documented. For documentation of the various functions, see the specific 
modules.                                                                    
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Import of modules                                                           #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

import import_files_module as impfil
import process_split_vep_module as provep
import process_variants_module as provar

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Defining the main function                                                  #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

def main():
    """
    This is the entry point of the Python program, and the execution of the
    program consists of calling this function (see line 208).
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
    
    ### Part 0.3: Get the list of patient-specific cram file locations from
    # file
    cram_file_list = impfil.read_cram_file_list(
            path_to_cram_list_file = str('file.txt'))
    
    ### Part 0.4: Get the appropriate cram file location from the patient ID
    path_to_patient_cram_file = impfil.get_path_to_cram_file_from_patient_ID(
            patient_ID = patient_ID,
            cram_file_list = cram_file_list)
    
    ### Part 0.5: Path to patient-specific folder and files
    # Assign path to patient-specific folder  and IGV screenshots file
        
    path_to_report_folder = str('path/on/server/')
    
    # Assign path to IGV batch file
    path_to_IGV_batch_file = str(path_to_report_folder + 
                                 'igv_snapshots/batch_file_' +
                                 str(patient_ID) +
                                 '.bat')
            
    # Assign path to snapshot directory
    path_to_snapshot_directory = str(path_to_report_folder + 
                                 'igv_snapshots/')
            
    ### Part 0.6: Assign path to VCF data
    
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
    # Part 2: Transformation and processing of the variant call dataframe   #
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
                    
    # Step 3.1: Processing of expanded gene list variants (groupID = group 2-6)
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    # Extraction of variants from group 2-6 (expanded gene list)
    # and quality filtering based on genotype quality (all variants with 
    # genotype quality < 50 are not considered)
    expanded_variants = provar.extract_group_list_variants(
                variant_calls_dataframe = variant_calls,
                groupID = 'grp2-6')
            
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Part 4: Generate IGV read mapping batch file                            #
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    # Generate batch file to be used in IGV to visualise reads
    provar.generate_IGV_batch_file(
            variant_calls_curated = curated_variants, 
            variant_calls_expanded = expanded_variants,
            path_to_IGV_batch_file = path_to_IGV_batch_file,
            path_to_snapshot_directory = path_to_snapshot_directory,
            path_to_cram_file = path_to_patient_cram_file)
    
    return 0

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Execution of the main function                                              #
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

if __name__=="__main__":
    main()

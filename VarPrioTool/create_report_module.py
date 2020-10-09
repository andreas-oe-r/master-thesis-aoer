"""
# -*- coding: utf-8 -*-

### File name:                                                              
create_report_module.py                              

### General description of the script:                                      
This modules contains the functions for creation of the PDF report.                                                                  
                                                                          
"""

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Imports                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# Description:
# The module Reportlab is used for creation of the PDF. The use of this module 
# is based on the "canvas" and "doc" objects. The PDF is created by adding
# "flowables" to the doc object.

from reportlab.platypus import (Flowable, Paragraph, SimpleDocTemplate, Spacer, 
                                PageBreak, Image)
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.lib.pagesizes import letter
from reportlab.platypus import Table, TableStyle
from reportlab.lib import colors
from reportlab.pdfbase import pdfmetrics
from reportlab.lib.units import inch

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Classes                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

class insert_horizontal_line(Flowable):
    """
    This function creates a horizontal line as a type Reportlab flowable.
        
    """
    #----------------------------------------------------------------------
    def __init__(self, width, height=0):
        Flowable.__init__(self)
        self.width = width
        self.height = height
    #----------------------------------------------------------------------
    def __repr__(self):
        return "Line(w=%s)" % self.width
    #----------------------------------------------------------------------
    def draw(self):
        """
        draw the line
        """
        self.canv.line(0, self.height, self.width, self.height)

# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
# Functions                                                                  
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

# Function 0: ...

def generate_pdf_object(path_to_pdf_file_str):
    """
    Generates the PDF object.

    Parameters
    ----------
    path_to_pdf_file_str : str
        The file location of the output PDF report
    """
    pdf_object = SimpleDocTemplate(path_to_pdf_file_str,
                            pagesize = letter,
                            rightMargin = 20,
                            leftMargin = 20,
                            topMargin = 20,
                            bottomMargin = 20)
    
    return pdf_object

# Function 1: ...

def footer(canvas, doc):
    """
    Generates the PDF object.

    Parameters
    ----------
    canvas : Reportlab canvas object
        The canvas on which the footer will be placed
        
    doc : Reportlab doc object

    """
    styles = getSampleStyleSheet()
    styleN = styles['Normal']
    canvas.saveState()
    P = Paragraph('<a href="#ANCHOR_1" color="blue">Return to top</a>',
                  styleN)
    w, h = P.wrap(doc.width, doc.bottomMargin)
    P.drawOn(canvas, doc.leftMargin, h)
    canvas.restoreState()


# Function 2: ...

def insert_summary_table(data_table, list_of_anchors):
    """
    This function creates the summary table.

    Parameters
    ----------
    data_table : list of lists
        The list of list containiing the data to be inserted
        
    list_of_anchors : list
        A list of hyperlinks
    
    Returns
    -------
    report_table : reportlab flowable object
        The reportlab flowable which must be inserted into the list of 
        flowables
    """
    # Replace the placeholder variant number with the link anchor
    for variant_row in data_table:
        variant_row[0] = list_of_anchors.pop(0)
    
    report_table = Table(data_table, hAlign = 'LEFT')
   
    report_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))

    report_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,-1),colors.whitesmoke)]))    
    report_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(-1,0),colors.lightgrey)]))

    return report_table


def insert_statistics_table(data_table):
    """
    This function creates the statistics table.

    Parameters
    ----------
    data_table : list of lists
        The list of list containiing the data to be inserted
        
    Returns
    -------
    report_table : reportlab flowable object
        The reportlab flowable which must be inserted into the list of 
        flowables
    """
    report_table = Table(data_table, hAlign = 'LEFT')
   
    report_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))

    report_table.setStyle(TableStyle([
            ('BACKGROUND',(1,0),(1,8),colors.whitesmoke)]))    
    report_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(0,8),colors.lightgrey)]))

    return report_table

def insert_expanded_info_table(data_table, anchor):
    """
    This function creates individual variant information boxes. Each box
    is made up of a number of individual tables. 

    Parameters
    ----------
    data_table : list of lists
        The list of list containiing the data to be inserted
        
    anchor : str
        A hyperlink from the summary table  
    
    Returns
    -------
    list_of_individual_tables : list of reportlab flowable objects
        Each object in the list of flowable objects must be added to the over-
        all list of flowable objects
    """
    list_of_individual_tables = []
        
    # Create table 0 (Variant header) and set variables
    table_0 = []
    table_0_str = data_table[0]
    table_0_str.append(anchor)
    table_0.append(table_0_str)
    table_0_table = Table(table_0, colWidths=[570], hAlign = 'LEFT')
    table_0_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.white)]))
    table_0_table.setStyle(TableStyle(
            [('FONT', (0, 0), (-1, -1), 'Times-Roman', 18)]))
    table_0_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_0_table)

    # Create table 1 (General variant information header) and set variables
    table_1 = []
    table_1_str = data_table[1]
    table_1.append(table_1_str)
    table_1_table = Table(table_1, colWidths=[570], hAlign = 'LEFT')
    table_1_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.pink)]))
    table_1_table.setStyle(TableStyle(
            [('FONT', (0, 0), (-1, -1), 'Times-Bold', 14)]))
    table_1_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_1_table)

    # Create table 2 (Chr:Pos | Gene symbol | Type | Genotype |
    # Worst consequence | GQ | DP | AD) and set variables
    table_2 = data_table[2]
    table_2_table = Table(table_2, colWidths=[90, 70, 80, 60, 130, 35, 35, 70],
                          hAlign = 'LEFT')
    table_2_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    table_2_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey)]))
    table_2_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_2_table)

    # Create table 3 (Feature type | Biotype | VARIANT_CLASS | 
    # gnomAD AF | Transcripts | Existing_variation) and set variables
    table_3 = data_table[3]
    table_3_table = Table(table_3, colWidths=[130, 130, 130, 180], 
                          hAlign = 'LEFT')
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey)]))
    table_3_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    table_3_table.setStyle(TableStyle(
            [('LINEABOVE', (0, 0), (-1, 0), 1.5, colors.black)]))
    list_of_individual_tables.append(table_3_table)
    
    # Create table 4 (CADD Phred | SIFT | PolyPhen | ClinVar significance) 
    # and set variables
    table_4 = data_table[4]
    table_4_table = Table(table_4, colWidths=[130, 130, 130, 180], 
                          hAlign = 'LEFT')
    table_4_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    table_4_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey)]))
    table_4_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    table_4_table.setStyle(TableStyle(
            [('LINEABOVE', (0, 0), (-1, 0), 1.5, colors.black)]))
    list_of_individual_tables.append(table_4_table)

    # Create table 3 (Feature type | Biotype | VARIANT_CLASS | 
    # gnomAD AF | Transcripts | Existing_variation) and set variables
    table_3 = data_table[5]
    table_3_table = Table(table_3, colWidths=[130, 130, 130, 180], 
                          hAlign = 'LEFT')
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey)]))
    table_3_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    table_3_table.setStyle(TableStyle(
            [('LINEABOVE', (0, 0), (-1, 0), 1.5, colors.black)]))
    list_of_individual_tables.append(table_3_table)

    # Create table 3 (Feature type | Biotype | VARIANT_CLASS | 
    # gnomAD AF | Transcripts | Existing_variation) and set variables
    table_3 = data_table[6]
    table_3_table = Table(table_3, colWidths=[130, 130, 130, 180], 
                          hAlign = 'LEFT')
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.whitesmoke)]))
    table_3_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, 0), colors.lightgrey)]))
    table_3_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    table_3_table.setStyle(TableStyle(
            [('LINEABOVE', (0, 0), (-1, 0), 1.5, colors.black)]))
    list_of_individual_tables.append(table_3_table)


    # Create table 5 (Sequence information) and set variables
    table_5 = []
    table_5_str = data_table[7]
    table_5.append(table_5_str)
    table_5_table = Table(table_5, colWidths=[570], hAlign = 'LEFT')
    table_5_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.pink)]))
    table_5_table.setStyle(TableStyle(
            [('FONT', (0, 0), (-1, -1), 'Times-Bold', 14)]))
    table_5_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_5_table)

    # Create table 6  and set variables
    table_6 = data_table[8]
    table_6_table = Table(table_6, colWidths=[190, 380], hAlign = 'LEFT')
    table_6_table.setStyle(TableStyle([
            ('BACKGROUND',(1,0),(1,8),colors.whitesmoke)]))    
    table_6_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(0,8),colors.lightgrey)]))
    table_6_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_6_table)

    # Create table 7 (Annotation information) and set variables
    table_7 = []
    table_7_str = data_table[9]
    table_7.append(table_7_str)
    table_7_table = Table(table_7, colWidths=[570], hAlign = 'LEFT')
    table_7_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.pink)]))
    table_7_table.setStyle(TableStyle(
            [('FONT', (0, 0), (-1, -1), 'Times-Bold', 14)]))
    table_7_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_7_table)

    # Create table 8
    table_8 = data_table[10]
    table_8_table = Table(table_8, colWidths=[190, 380], hAlign = 'LEFT')
    table_8_table.setStyle(TableStyle([
            ('BACKGROUND',(1,0),(1,9),colors.whitesmoke)]))    
    table_8_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(0,9),colors.lightgrey)]))
    table_8_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_8_table)

    # Create table 9 (Annotation information)
    table_9 = []
    table_9_str = data_table[11]
    table_9.append(table_9_str)
    table_9_table = Table(table_9, colWidths=[570], hAlign = 'LEFT')
    table_9_table.setStyle(TableStyle(
            [('BACKGROUND', (0, 0), (-1, -1), colors.pink)]))
    table_9_table.setStyle(TableStyle(
            [('FONT', (0, 0), (-1, -1), 'Times-Bold', 14)]))
    table_9_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_9_table)

    # Create table 10
    table_10 = data_table[12]
    table_10_table = Table(table_10, colWidths=[190, 380], hAlign = 'LEFT')
    table_10_table.setStyle(TableStyle([
            ('BACKGROUND',(1,0),(1,9),colors.whitesmoke)]))    
    table_10_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(0,9),colors.lightgrey)]))
    table_10_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    list_of_individual_tables.append(table_10_table)

    # Create table 11
    table_11 = data_table[13]
    table_11_table = Table(table_11, colWidths=[190, 380], hAlign = 'LEFT')
    table_11_table.setStyle(TableStyle([
            ('BACKGROUND',(1,0),(1,9),colors.whitesmoke)]))    
    table_11_table.setStyle(TableStyle([
            ('BACKGROUND',(0,0),(0,9),colors.lightgrey)]))
    table_11_table.setStyle(TableStyle(
            [("BOX", (0, 0), (-1, -1), 0.25, colors.black),
            ('INNERGRID', (0, 0), (-1, -1), 0.25, colors.black)]))
    table_11_table.setStyle(TableStyle(
            [('LINEBELOW', (0, 0), (-1, 0), 1.5, colors.black)]))
    list_of_individual_tables.append(table_11_table)

    return list_of_individual_tables

# Function 2: 

def generate_variant_analysis_report(patient_id, pdf_object, pdf_elements,
                                     path_to_IGV_snapshots_directory):
    """
    This function creates the actual PDF, using the other functions and various
    input. 

    Parameters
    ----------
    patient_id : str
        The patient ID to be used in the header
        
    pdf_object : Reportlab doc object
        The reportlab doc object that has been initialised  
        
    pdf_elements : List of pandas dataframes
        All of the data that is to be inserted into the pdf report
        
    path_to_IGV_snapshots_directory : str
        The path to the IGV snapshots directroy
    
    """
    
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Set various paragraph styles and create list for PDF contents
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    styles = getSampleStyleSheet()
    
    # Set main header style
    styles.add(ParagraphStyle(name='main_header', alignment=TA_JUSTIFY))
    styles["main_header"].fontSize = 20 
    styles["main_header"].leading = 14

    # Set sub-header 1 style
    styles.add(ParagraphStyle(name='sub_header_1', alignment=TA_JUSTIFY))
    styles["sub_header_1"].fontSize = 16 
    styles["sub_header_1"].leading = 14

    # Set sub-header 2 style
    styles.add(ParagraphStyle(name='sub_header_2', alignment=TA_JUSTIFY))
    styles["sub_header_2"].fontSize = 16 
    styles["sub_header_2"].leading = 14

    # Set sub-header 3 style
    styles.add(ParagraphStyle(name='sub_header_3', alignment=TA_JUSTIFY))
    styles["sub_header_3"].fontSize = 14 
    styles["sub_header_3"].leading = 14
    
    # Set sub-header 4 style
    styles.add(ParagraphStyle(name='sub_header_4', alignment=TA_JUSTIFY))
    styles["sub_header_4"].fontSize = 14 
    styles["sub_header_4"].leading = 14
    styles["sub_header_4"].bulletIndent = 1
    
    # Set normal style
    styles.add(ParagraphStyle(name='normal', alignment=TA_JUSTIFY))
    styles["normal"].fontSize = 10 
    styles["normal"].leading = 10

    flowables = []
    
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert main header and part 1 (summary of results) header
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    # Insert line before main header
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))

    # Add main header
    text = '''<a name="ANCHOR_1"/><font color="black">Thromb variant 
    analysis report for {}</font>'''.format(patient_id)
    para = Paragraph(text, style = styles["main_header"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
            
    # Add sub header 1 beneath main header
    text = '''Results from gene list based variant prioritization relevant for
            probands'''
    para = Paragraph(text, style = styles["sub_header_1"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 6))

    text = '''who exhibit an abnormal bleeding phenotype'''
    para = Paragraph(text, style = styles["sub_header_1"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 18))

    # Insert line after main header
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))

    # Add header for part 1
    text = "Part 1/2 (summary of results)"
    para = Paragraph(text, style = styles["sub_header_2"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))

    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert part 1 (summary of results) for curated gene list
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    # Add header for part 1, curated gene list
    text = "1. Curated gene list / gene panel (diagnostic-grade genes)"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
    
    # Insert line after subheader
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))

    # Add header for part 1 – 1.1.
    text = "1.1. Summary of variants found"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))

    # Code for including hyperlinks in the table (needs complete rework of 
    # genereate_summary_table and insert_summary table, such that a list of
    # hyperlink anchors can be passed to generate_summary_table)
    
    # Create a list that will hold the variant ID's as a paragraph object
    # thereby enabling hyperlinking
    list_of_variant_anchors = []
    list_of_variant_anchors.append('VarNum')
    
    # Iterate through the summary list of variants (first row is column head-
    # ings), so the number of variants is the length minus 1
    for i in range(len(pdf_elements[0]) - 1):
        var_text = str('Var_' + str(i))
        anchor_text = str('ANCHOR_CUR_SUM_' + str(i))
        text = str(
                '<a href="#' + anchor_text + '" color="blue">' + var_text +
                '</a>')
        L = pdfmetrics.stringWidth('approx_len', "Times-Roman", 10)
        inside_Table = Table(
                [[Paragraph(text, style = styles["normal"])]], colWidths=L)
        list_of_variant_anchors.append(inside_Table)

    # Insert the summary table
    report_table = insert_summary_table(pdf_elements[0], list_of_variant_anchors)
    flowables.append(report_table)
    flowables.append(Spacer(1, 14))

    # Add header for part 1 – 1.2.
    text = "1.2. General statistics"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))

    # Insert the statistics table for curated gene list variants
    report_table = insert_statistics_table(pdf_elements[1])
    flowables.append(report_table)
    flowables.append(Spacer(1, 14))

    # Add page break
    flowables.append(PageBreak())

    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert part 1 (summary of results) for expanded gene list
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    # Add header for part 1, expanded gene list
    text = "2. Expanded gene list"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
    
    # Insert line after subheader
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))

    # Add header for part 1 – 2.1.
    text = "2.1. Summary of variants found"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))

    # Create a list that will hold the variant ID's as a paragraph object
    # thereby enabling hyperlinking
    list_of_variant_anchors = []
    list_of_variant_anchors.append('VarNum')
    
    # Iterate through the summary list of variants (first row is column head-
    # ings), so the number of variants is the length minus 1
    for i in range(len(pdf_elements[2]) - 1):
        var_text = str('Var_' + str(i))
        anchor_text = str('ANCHOR_EXP_SUM_' + str(i))
        text = str(
                '<a href="#' + anchor_text + '" color="blue">' + var_text +
                '</a>')
        L = pdfmetrics.stringWidth('approx_len', "Times-Roman", 10)
        inside_Table = Table(
                [[Paragraph(text, style = styles["normal"])]], colWidths=L)
        list_of_variant_anchors.append(inside_Table)

    # Insert the summary table
    report_table = insert_summary_table(pdf_elements[2], list_of_variant_anchors)
    flowables.append(report_table)
    flowables.append(Spacer(1, 14))

    # Add header for part 1 – 2.2.
    text = "2.2. General statistics"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
   
    # Insert the statistics table for expanded gene list variants
    report_table = insert_statistics_table(pdf_elements[3])
    flowables.append(report_table)
    flowables.append(Spacer(1, 14))

    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert part 2 header
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #

    # Add page break
    flowables.append(PageBreak())
    
    # Add header for part 2 
    text = "Part 2/2 (expanded information)"
    para = Paragraph(text, style = styles["sub_header_2"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))

    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert part 2 (variant information) for curated gene list
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    # Add header for part 2 – 3.
    text = "3. Curated gene list / gene panel (diagnostic-grade genes)"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
    
    # Insert line after subheader
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))
        
    # Insert table for each variant
    total_num = 0
    num = 0
    for i in pdf_elements[4]:
        if num % 1 == 0 and num != 0:
            flowables.append(PageBreak())
        anchor_text = str('ANCHOR_CUR_SUM_' + str(num))
        var_text = str('Variant number ' + str(num))
        text = str('<a name="' + anchor_text + '"/><font color="black">' +
                   var_text + ' (curated gene list) </font>')
        L = pdfmetrics.stringWidth('approximate_length_is_this_long_and_a_bit_more', "Times-Bold", 14)
        inside_Table = Table(
                [[Paragraph(text, style = styles["sub_header_4"])]], colWidths=L)
        report_tables = insert_expanded_info_table(i, inside_Table)
        for i in report_tables:
            flowables.append(i)
        flowables.append(Spacer(1, 14))
        flowables.append(Spacer(1, 14))
        snapshot = Image(str(path_to_IGV_snapshots_directory + str(total_num) + '.png'),
                         7*inch, 7*inch)
        flowables.append(snapshot)
        flowables.append(Spacer(1, 14))
        num = num + 1
        total_num = total_num + 1
    
    # Add page break
    flowables.append(PageBreak())
    
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    # Insert part 2 (variant information) for expanded gene list
    # ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––– #
    
    # Add header for part 2 – 4.
    text = "4. Expanded gene list"
    para = Paragraph(text, style = styles["sub_header_3"]) 
    flowables.append(para)
    flowables.append(Spacer(1, 14))
    
    # Insert line after subheader
    line = insert_horizontal_line(500)
    flowables.append(line)
    flowables.append(Spacer(1, 14))

    # Insert table for each variant
    num = 0
    for i in pdf_elements[5]:
        if num % 1 == 0 and num != 0:
            flowables.append(PageBreak())
        anchor_text = str('ANCHOR_EXP_SUM_' + str(num))
        var_text = str('Variant number ' + str(num))
        text = str('<a name="' + anchor_text + '"/><font color="black">' +
                   var_text + ' (expanded gene list) </font>')
        L = pdfmetrics.stringWidth('approximate_length_is_this_long_and_a_bit_more_so_actually_this_long', "Times-Bold", 14)
        inside_Table = Table(
                [[Paragraph(text, style = styles["sub_header_4"])]], colWidths=L)
        report_tables = insert_expanded_info_table(i, inside_Table)
        for i in report_tables:
            flowables.append(i)
        flowables.append(Spacer(1, 14))
        flowables.append(Spacer(1, 14))
        snapshot = Image(str(path_to_IGV_snapshots_directory + str(total_num) + '.png'),
                         7*inch, 7*inch)
        flowables.append(snapshot)
        flowables.append(Spacer(1, 14))
        num = num + 1
        total_num = total_num + 1

    pdf_object.build(flowables,
                     onLaterPages = footer)

    return 0

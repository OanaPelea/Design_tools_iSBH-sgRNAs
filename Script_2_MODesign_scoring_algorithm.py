#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:08:24 2020

@author: oanapelea
"""
from iSBH_functions import determine_mfe_structure
from iSBH_functions import return_percentage_similarity_2_sequences
from iSBH_functions import get_DNA_reverse_complement
from iSBH_functions import openess_mfe_structure
from iSBH_functions import standard_sequence_format

input_filename='intermediate_output/MODesign_sgRNA_1_Alice_eRNA_trigger_25_05_2023.txt'

output_file=open('outputs/1_sorted_MODesign_sgRNA_1_Alice_eRNA_trigger_25_05_2023.csv', 'w')

original_spacer_sequence='TGGGAAGCAGAGGCAAAGGG'

trigger_name='Alice_eRNA'
full_trigger_sequence=standard_sequence_format('CTTTGACTCCTCAAACAGAATTAAGCATTGTCCATCTCCCCTTCCAGAGACTGGGTCTTCTTTGTATCCCTCACACTCAGCTGAGGGCCTGGAGCAATGGAGATACTCAAAAAGCCCTCAAGAATGAAGGGAGGGATGAATATTATATGAGATAAGAAGCAAGAAAACGTGATTCCTGACACA')

trigger_secondary_structure=determine_mfe_structure(full_trigger_sequence,"rna")
print trigger_secondary_structure


for line in open(input_filename):
    f=line.split()
    design_name=f[0][0:19]
    iSBH_sequence=f[1]
    truncated_spacer_sequence=f[2]
    folding_probability=float(f[3])
    long_trigger_sub_sequence=''
    small_trigger_sub_sequence=''
    
    if design_name=='Sl14_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:44])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
    
    if design_name=='Sl16_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:46])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl16.2_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:(46+5)])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl18_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:48])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl18.2_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:(48+5)])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20] 
      
    if design_name=='Sl20_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:50])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl22_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:52])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl22.2_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:(52+5)])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl24_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:54])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl24.2_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:(54+5)])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl26_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:56])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl28_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:58])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]

    if design_name=='Sl30_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:60])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
      
    if design_name=='Sl32_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:62])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20] 
      
    if design_name=='Sl34_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:64])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]      
      
    if design_name=='Sl36_giraffe_design':
      long_trigger_sub_sequence=get_DNA_reverse_complement(iSBH_sequence[0:66])
      small_trigger_sub_sequence=long_trigger_sub_sequence[0:-20]
          
    if small_trigger_sub_sequence in full_trigger_sequence:
       trigger_start_index=full_trigger_sequence.index(small_trigger_sub_sequence)
       small_trigger_end_index=trigger_start_index+len(small_trigger_sub_sequence)
       long_trigger_end_index=trigger_start_index+len(long_trigger_sub_sequence)
      
      
       extended_trigger_sequence_including_spacer_star=full_trigger_sequence[trigger_start_index:long_trigger_end_index]   

       percentage_trigger_similarity=return_percentage_similarity_2_sequences(extended_trigger_sequence_including_spacer_star, long_trigger_sub_sequence)
       percentage_spacer_similarity=return_percentage_similarity_2_sequences(truncated_spacer_sequence, original_spacer_sequence)
       secondary_structure_small_trigger_sub_sequence=trigger_secondary_structure[trigger_start_index:small_trigger_end_index]
       openess_secondary_structure_small_trigger_sub_sequence=openess_mfe_structure(secondary_structure_small_trigger_sub_sequence)
      
       if percentage_trigger_similarity!='error':
        
         quality_score=folding_probability*folding_probability*folding_probability*percentage_spacer_similarity*float(percentage_trigger_similarity)*openess_secondary_structure_small_trigger_sub_sequence
         output=trigger_name+','+design_name+','+iSBH_sequence+','+str(quality_score)+','+small_trigger_sub_sequence+',,,\n'
         output_file.write(output)
output_file.close()
#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 30 09:26:05 2020

@author: oanapelea
"""

from iSBH_functions import generate_symbolic_iSBH_fold_from_given_iSBH_components_no_toehold
from iSBH_functions import generate_iSBH_symbolic_fold_from_given_inputs_with_toehold
from iSBH_functions import generate_iSBH_sensing_modules_for_given_trigger_and_specifications
from iSBH_functions import standard_sequence_format
from iSBH_functions import determine_most_stable_backfold_for_given_spacer
from iSBH_functions import identifying_corresponding_MM_positions
from iSBH_functions import assemble_preliminary_iSBH_sequence_no_bulges
from iSBH_functions import asemble_final_iSBHs_from_iSBHs_without_bulges_and_variable_bulge_sizes
from iSBH_functions import generate_iSBH_sensing_modules_for_given_trigger_and_specifications
from iSBH_functions import append_toeholds_to_iSBHs_with_bulges
from iSBH_functions import determine_consecutive_base_pairing_interactions_between_iSBHs_no_toehold_and_trigger
from iSBH_functions import selecting_no_toehold_iSBHs_that_still_maintain_final_folding_when_adding_tracr
from iSBH_functions import get_DNA_reverse_complement
from iSBH_functions import return_all_bulge_combinations_except_input
from iSBH_functions import determine_mfe_structure
from iSBH_functions import determine_probability_secondary_structure
from iSBH_functions import return_percentage_similarity_2_sequences
from iSBH_functions import standard_sequence_format


#input several features for the iSBHs you would like to design
folding_probability_cut_off=0.4
today_date='25_05_2023'

list_design_components=[['..(((..(((', ')))..)))..', 14,'Sl14_giraffe_design'], ['..(((..(((', ')))..)))..', 16,'Sl16_giraffe_design'], ['..(((..(((', ')))..)))..', 18,'Sl18_giraffe_design'], ['..(((..(((', ')))..)))..', 20,'Sl20_giraffe_design'],['..(((..(((', ')))..)))..', 22,'Sl22_giraffe_design'], ['..(((..(((', ')))..)))..', 24,'Sl24_giraffe_design'], ['..(((..(((', ')))..)))..', 26,'Sl26_giraffe_design'], ['..(((..(((', ')))..)))..', 28,'Sl28_giraffe_design'], ['..(((..(((', ')))..)))..', 30,'Sl30_giraffe_design'], ['..(((..(((', ')))..)))..', 32,'Sl32_giraffe_design']]

#lists of spacers and triggers for which I should carry out the original testing
list_of_trigger_names=['Alice_eRNA']

list_of_triggers=[standard_sequence_format('CTTTGACTCCTCAAACAGAATTAAGCATTGTCCATCTCCCCTTCCAGAGACTGGGTCTTCTTTGTATCCCTCACACTCAGCTGAGGGCCTGGAGCAATGGAGATACTCAAAAAGCCCTCAAGAATGAAGGGAGGGATGAATATTATATGAGATAAGAAGCAAGAAAACGTGATTCCTGACACA')]
input_filename='outputs/optimal_spacer_folds.csv'

spacer_star_fold='((((((((((..(((..((('
spacer_fold=')))..)))..))))))))))'

for line in open(input_filename).readlines():
  if 'sgRNA_name' not in line:
    f=line.split(',')
    sgRNA_name=f[0]
    spacer_sequence=f[1]
    spacer_star_sequence=f[2] 
    original_15_nt_spacer_star_sequece=spacer_star_sequence[0:15]
    spacer_sequence_17nt=spacer_sequence[3:]
    
    print original_15_nt_spacer_star_sequece, spacer_sequence_17nt
  
    index=0
    for trigger in list_of_triggers:
        output_filename='intermediate_output/'+'MODesign_'+sgRNA_name+'_'+list_of_trigger_names[index]+"_trigger_"+today_date+'.txt'
        output_file=open(output_filename,'w')
        print list_of_trigger_names[index]
        index=index+1 
    
        for design_specifications in list_design_components:
            giraffe_star_fold= design_specifications[0]
            giraffe_fold=design_specifications[1]
            loop_length=design_specifications[2]
            design_name=design_specifications[3]

            desired_fold=generate_symbolic_iSBH_fold_from_given_iSBH_components_no_toehold(giraffe_star_fold, giraffe_fold,loop_length)[0]
            desired_fold_iSBH_with_bulges_only_in_the_spacer=spacer_star_fold+('('*len(giraffe_star_fold))+("."*loop_length)+(')'*len(giraffe_fold))+spacer_fold
    
            fully_sensing_trigger_length=int(generate_symbolic_iSBH_fold_from_given_iSBH_components_no_toehold(giraffe_star_fold, giraffe_fold,loop_length)[1])
            trigger_length=20+fully_sensing_trigger_length
     
            list_sensing_modules=generate_iSBH_sensing_modules_for_given_trigger_and_specifications(trigger, trigger_length)

            list_iSBHs_no_bulges=[]
     
            #generating a list of native sequences that would form the loop in a desired way
            for sensing_module in list_sensing_modules:
                spacer_and_giraffe_extension=get_DNA_reverse_complement(sensing_module[:-loop_length])
                initial_iSBH_endogenous_sequence=sensing_module+spacer_and_giraffe_extension
        

                modular_iSBH=original_15_nt_spacer_star_sequece+initial_iSBH_endogenous_sequence[15:-17]+spacer_sequence_17nt
 
                mfe_structure_modular_iSBH=determine_mfe_structure(modular_iSBH, 'rna')
                
        
                if mfe_structure_modular_iSBH == desired_fold_iSBH_with_bulges_only_in_the_spacer:
                  list_iSBHs_no_bulges.append(modular_iSBH)
            
            if len(list_iSBHs_no_bulges) > 0:
               output_designing_iSBHs_with_bulges=asemble_final_iSBHs_from_iSBHs_without_bulges_and_variable_bulge_sizes(list_iSBHs_no_bulges,desired_fold)
               if output_designing_iSBHs_with_bulges!= 'NA':
          
                  list_iSBH_sequences_with_bulges=output_designing_iSBHs_with_bulges[0]
                  list_folding_probabilities_iSBH_sequences_with_bulges=output_designing_iSBHs_with_bulges[1]
            
                  #outputting iSBH characteristics
                 
                  index_2=0
                  list_unique_iSBH_designs=[]
                  while index_2 < len(list_folding_probabilities_iSBH_sequences_with_bulges):
                                probability=float(list_folding_probabilities_iSBH_sequences_with_bulges[index_2])
                                if probability > folding_probability_cut_off:
                                      iSBH_sequence=list_iSBH_sequences_with_bulges[index_2]
                                      info_trigger_sub_sequence=determine_consecutive_base_pairing_interactions_between_iSBHs_no_toehold_and_trigger(iSBH_sequence, trigger, trigger_length)
                                      length_trigger_sub_sequence=info_trigger_sub_sequence[1]
                                      trigger_sub_sequence=info_trigger_sub_sequence[2]
                                      folding_probability=str(probability)
                                      
                                      if iSBH_sequence not in list_unique_iSBH_designs:
                                        list_unique_iSBH_designs.append(iSBH_sequence)
                                        output=str(design_name)+" "+str(iSBH_sequence)+" "+str(iSBH_sequence[-20:])+" "+str(folding_probability)+" "+ str(length_trigger_sub_sequence)+' '+str(trigger_sub_sequence)+"   \n"
                                        output_file.write(output)
                              
                                        print output
                                      
                                    
                                index_2=index_2+1

                
        output_file.close()
    
     
     

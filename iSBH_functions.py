#!/usr/bin/bash python2
# -*- coding: utf-8 -*-
"""
Created on Tue Oct  9 22:30:51 2018

@author: oanapelea
"""
import random
import os
import stat
import subprocess
import numpy as np

subprocess.call('export NUPACKHOME=/Users/oanapelea/Documents/useful_programs/nupack3.2.2', shell=True)

def separate_iSBH_components(name, sequence):
    if "spa" in name:
        spacer_star=''
        hairpin_star=''
        loop=''
        hairpin=''
        spacer=sequence
    if (("SL14" in name) and ("original_design" in name)):
        if len(sequence) == 54:
          spacer_star=sequence[0:20]
          hairpin_star=''
          loop=sequence[20:34]
          hairpin=''
          spacer=sequence[34:54]
        if len(sequence) == 52:           
          spacer_star=sequence[0:19]
          hairpin_star=''
          loop=sequence[19:33]
          hairpin=''
          spacer=sequence[33:52]
    if (("SL20" in name) and ("original_design" in name)):
        if len(sequence) == 60:         
          spacer_star=sequence[0:20]
          hairpin_star=''
          loop=sequence[20:40]
          hairpin=''
          spacer=sequence[40:60]
        if len(sequence) == 58:  
          spacer_star=sequence[0:19]
          hairpin_star=''
          loop=sequence[19:39]
          hairpin=''
          spacer=sequence[39:58]   
    if (("SL25" in name) and ("original_design" in name)):
        if len(sequence) == 65: 
          spacer_star=sequence[0:20]
          hairpin_star=''
          loop=sequence[20:45]
          hairpin=''
          spacer=sequence[45:65]
        if len(sequence) == 63: 
          spacer_star=sequence[0:19]
          hairpin_star=''
          loop=sequence[19:44]
          hairpin=''
          spacer=sequence[44:63]
    if ("SL14" in name) and ("giraffe" in name):
        if len(sequence) == 74:
          spacer_star=sequence[0:20]
          hairpin_star=sequence[20:30]
          loop=sequence[30:44]
          hairpin=sequence[44:54]
          spacer=sequence[54:74]
        if len(sequence) == 72:           
          spacer_star=sequence[0:19]
          hairpin_star=sequence[19:29]
          loop=sequence[29:43]
          hairpin=sequence[43:53]
          spacer=sequence[53:72]
    if ("SL20" in name) and ("giraffe" in name):
        if len(sequence) == 80:         
          spacer_star=sequence[0:20]
          hairpin_star=sequence[20:30]
          loop=sequence[30:50]
          hairpin=sequence[50:60]
          spacer=sequence[60:80]
        if len(sequence) == 78:  
          spacer_star=sequence[0:19]
          hairpin_star=sequence[19:29]
          loop=sequence[19:49]
          hairpin=sequence[49:59]
          spacer=sequence[59:78]   
    if ("SL25" in name) and ("giraffe" in name):
        if len(sequence) == 85: 
          spacer_star=sequence[0:20]
          hairpin_star=sequence[20:30]
          loop=sequence[30:55]
          hairpin=sequence[55:65]
          spacer=sequence[65:85]
        if len(sequence) == 83: 
          spacer_star=sequence[0:19]
          hairpin_star=sequence[19:29]
          loop=sequence[29:54]
          hairpin=sequence[54:64]
          spacer=sequence[64:83]
    return (name, spacer_star, hairpin_star, loop, hairpin, spacer)
    
def create_random_DNA_sequence(desired_length):
    length=int(desired_length)
    nucleotide_list=['A','C','T','G']
    random_sequence=''
    index=0
    while index < length:
        random_nucleotide=random.choice(nucleotide_list)
        random_sequence=random_sequence+random_nucleotide
        index=index+1
    return (random_sequence)

def create_random_DNA_sequence_considering_GC(desired_length, GC_min, GC_max):
    length=float(desired_length)
    nucleotide_list=['A','C','T','G']
    sequences_list=[]
    while len(sequences_list) <=3:
      sequence=''
      index=0
      GC_content=0
      while index < length:
        random_nucleotide=random.choice(nucleotide_list)
        sequence=sequence+random_nucleotide
        if random_nucleotide == 'G' or random_nucleotide=='C':
          GC_content=GC_content+1
        index=index+1        
      GC_percentage=float((GC_content/length)*100)
      if ((GC_percentage >= float(GC_min)) and (GC_percentage <= float(GC_max))):
        sequences_list.append(sequence) 
    desired_sequence=sequences_list[1]
    return (desired_sequence)

def assess_GC_content(sequence):
    sequence_length=float(len(sequence))
    index=0
    GC_content=0
    while index<sequence_length:
        nucleotide=sequence[index]
        if nucleotide == "G" or nucleotide=='C':
            GC_content=GC_content+1
        index=index+1
    GC_percentage=float((GC_content/sequence_length)*100)
    return (GC_percentage)

def get_DNA_reverse_complement(DNA_sequence):
    sequence_length=len(DNA_sequence)
    index=0
    complementary_reverse_sequence=''
    while index < sequence_length:
        nucleotide=DNA_sequence[index]
        if nucleotide=="A":
            complementary_nucleotide="T"
        if nucleotide=="C":
            complementary_nucleotide="G"
        if nucleotide=="T":
            complementary_nucleotide="A"
        if nucleotide=="G":
            complementary_nucleotide="C"
        complementary_reverse_sequence=complementary_nucleotide+complementary_reverse_sequence
        index=index+1
    return(complementary_reverse_sequence)
    
def get_RNA_reverse_complement(RNA_sequence):
    sequence_length=len(RNA_sequence)
    index=0
    complementary_reverse_sequence=''
    while index < sequence_length:
        nucleotide=RNA_sequence[index]
        if nucleotide=="A":
            complementary_nucleotide="U"
        if nucleotide=="C":
            complementary_nucleotide="G"
        if nucleotide=="U":
            complementary_nucleotide="A"
        if nucleotide=="G":
            complementary_nucleotide="C"
        complementary_reverse_sequence=complementary_nucleotide+complementary_reverse_sequence
        index=index+1
    return(complementary_reverse_sequence)

def get_mismatched_RNA_sequence(sequence,MM_position_1, MM_position_2, MM_position_3, MM_position_4):
    sequence_length=len(sequence)
    MM_list=[]
    if (MM_position_1 !=0) and (int(MM_position_1-1)<=sequence_length):
       MM_list.append(int(MM_position_1-1))
    if (MM_position_2 !=0) and (int(MM_position_2-1)<=sequence_length):
       MM_list.append(int(MM_position_2-1)) 
    if (MM_position_3 !=0) and (int(MM_position_3-1)<=sequence_length):
       MM_list.append(int(MM_position_3-1))
    if (MM_position_4 !=0) and (int(MM_position_4-1)<=sequence_length):
       MM_list.append(int(MM_position_4-1)) 
    index=0
    mismatched_sequence=''
    while index < sequence_length:
      nucleotide=sequence[index]
      if index in MM_list:
          if nucleotide=='A':
              changed_nucleotide=random.choice(["U","C", "G"])
          if nucleotide=='U':
              changed_nucleotide=random.choice(["G","A","C"])
          if nucleotide=='G':
              changed_nucleotide=random.choice(["A","U","C"])
          if nucleotide=='C':
              changed_nucleotide=random.choice(["G","A", "U"])
          mismatched_sequence=mismatched_sequence+changed_nucleotide
      else:
          mismatched_sequence=mismatched_sequence+nucleotide
      index=index+1
    
    return (mismatched_sequence)


def get_DNA_mismatched_reverse_complement(DNA_seq,MM_position_1, MM_position_2, MM_position_3, MM_position_4):    
  sequence_length=len(DNA_seq)
  MM_list=[]
  if (MM_position_1 !="NA") and (int(MM_position_1-1)<=sequence_length):
     MM_list.append(int(MM_position_1-1))
  if (MM_position_2 !="NA") and (int(MM_position_2-1)<=sequence_length):
     MM_list.append(int(MM_position_2-1)) 
  if (MM_position_3 !="NA") and (int(MM_position_3-1)<=sequence_length):
     MM_list.append(int(MM_position_3-1))
  if (MM_position_4 !="NA") and (int(MM_position_4-1)<=sequence_length):
     MM_list.append(int(MM_position_4-1)) 
  index=0
  mismatched_sequence=''
  while index < sequence_length:
      nucleotide=DNA_seq[index]
      if index in MM_list:
          if nucleotide=='A':
              changed_nucleotide=random.choice(["T","C", "G"])
          if nucleotide=='G':
              changed_nucleotide=random.choice(["A","T","C"])
          if nucleotide=='T':
              changed_nucleotide=random.choice(["G","A","C"])
          if nucleotide=='C':
              changed_nucleotide=random.choice(["G","A", "T"])
          mismatched_sequence=mismatched_sequence+changed_nucleotide
      else:
          mismatched_sequence=mismatched_sequence+nucleotide
      index=index+1
      
  mismatched_sequence_reverse_complement=''
  index=0
  while index < sequence_length:
      nucleotide=mismatched_sequence[index]
      if nucleotide=="A":
        complementary_nucleotide="T"
      if nucleotide=="C":
        complementary_nucleotide="G"
      if nucleotide=="T":
        complementary_nucleotide="A"
      if nucleotide=="G":
        complementary_nucleotide="C"
      mismatched_sequence_reverse_complement=complementary_nucleotide+mismatched_sequence_reverse_complement
      index=index+1
  return (mismatched_sequence_reverse_complement)
  
def get_RNA_mismatched_reverse_complement(RNA_seq,MM_position_1, MM_position_2, MM_position_3, MM_position_4):    
  sequence_length=len(RNA_seq)
  MM_list=[]
  if (MM_position_1 !="NA") and (int(MM_position_1-1)<sequence_length):
     MM_list.append(int(MM_position_1-1))
  if (MM_position_2 !="NA") and (int(MM_position_2-1)<sequence_length):
     MM_list.append(int(MM_position_2-1)) 
  if (MM_position_3 !="NA") and (int(MM_position_3-1)<sequence_length):
     MM_list.append(int(MM_position_3-1))
  if (MM_position_4 !="NA") and (int(MM_position_4-1)<sequence_length):
     MM_list.append(int(MM_position_4-1)) 
  index=0
  mismatched_sequence=''
  while index < sequence_length:
      nucleotide=RNA_seq[index]
      if index in MM_list:
          if nucleotide=='A':
              changed_nucleotide=random.choice(["U","C", "G"])
          if nucleotide=='G':
              changed_nucleotide=random.choice(["A","U","C"])
          if nucleotide=='U':
              changed_nucleotide=random.choice(["G","A","C"])
          if nucleotide=='C':
              changed_nucleotide=random.choice(["G","A", "T"])
          mismatched_sequence=mismatched_sequence+changed_nucleotide
      else:
          mismatched_sequence=mismatched_sequence+nucleotide
      index=index+1
      
  mismatched_sequence_reverse_complement=''
  index=0
  while index < sequence_length:
      nucleotide=mismatched_sequence[index]
      if nucleotide=="A":
        complementary_nucleotide="U"
      if nucleotide=="C":
        complementary_nucleotide="G"
      if nucleotide=="U":
        complementary_nucleotide="A"
      if nucleotide=="G":
        complementary_nucleotide="C"
      mismatched_sequence_reverse_complement=complementary_nucleotide+mismatched_sequence_reverse_complement
      index=index+1
  return (mismatched_sequence_reverse_complement)

def standard_sequence_format(any_sequence):
    sequence_length=len(any_sequence)
    standard_sequence_format=''
    index=0
    while index < sequence_length:
        nucleotide=any_sequence[index]
        if nucleotide=='a':
            standard_nucleotide='A'
        if nucleotide=='t':
            standard_nucleotide='T'
        if nucleotide=='u':
            standard_nucleotide='U'        
        if nucleotide=='g':
            standard_nucleotide='G'
        if nucleotide=='c':
            standard_nucleotide='C'
        if nucleotide=='A':
            standard_nucleotide='A'
        if nucleotide=='T':
            standard_nucleotide='T'
        if nucleotide=='U':
            standard_nucleotide='U'        
        if nucleotide=='G':
            standard_nucleotide='G'
        if nucleotide=='C':
            standard_nucleotide='C'

        standard_sequence_format=standard_sequence_format+standard_nucleotide
        index=index+1
    return(standard_sequence_format)

def reverse_sequence(any_sequence):
    sequence_length=len(any_sequence)
    index=0
    inverted_sequence=''
    while index < sequence_length:
        nucleotide=any_sequence[index]
        inverted_sequence=nucleotide+inverted_sequence
        index=index+1
    return(inverted_sequence)

def transcribe_sequence(DNA_sequence):
    RNA_sequence=''
    sequence_length=len(DNA_sequence)
    index=0
    while index < sequence_length:
        nucleotide=DNA_sequence[index]
        if nucleotide=='A':
            transcribed_nucleotide='A'
        if nucleotide=='T':
            transcribed_nucleotide='U'       
        if nucleotide=='G':
            transcribed_nucleotide='G'
        if nucleotide=='C':
            transcribed_nucleotide='C'
        RNA_sequence=RNA_sequence+transcribed_nucleotide
        index=index+1
    return(RNA_sequence)

def reverse_transcribe_sequence(RNA_sequence):
    DNA_sequence=''
    sequence_length=len(RNA_sequence)
    index=0
    while index < sequence_length:
        nucleotide=RNA_sequence[index]
        if nucleotide=='A':
            reversed_transcribed_nucleotide='A'
        if nucleotide=='U':
            reversed_transcribed_nucleotide='T'       
        if nucleotide=='G':
            reversed_transcribed_nucleotide='G'
        if nucleotide=='C':
            reversed_transcribed_nucleotide='C'
        DNA_sequence=DNA_sequence+reversed_transcribed_nucleotide
        index=index+1
    return(DNA_sequence)
    
def openess_mfe_structure(structure):
    unpaired_nucleotides=float(structure.count("."))
    total_nucleotides=float(len(structure))
    percentage_unpaired_nucleotides=float((unpaired_nucleotides/total_nucleotides)*100)
    return percentage_unpaired_nucleotides

def determine_mfe_structure(sequence, material): 
    sequence_length=len(sequence)
    out_file=open('intermediate_output/intermediate_output.in', "w")
    out_file.write(standard_sequence_format(sequence))
    out_file.write("\n")
    out_file.close()
    
    st = os.stat('intermediate_output/intermediate_output.in')
    os.chmod('intermediate_output/intermediate_output.in', st.st_mode | stat.S_IEXEC)
    
    inp_file='intermediate_output/intermediate_output'
    bash_command='mfe -material '+material+' ' + inp_file
    subprocess.call(bash_command, shell=True)
    
    for line in open('intermediate_output/intermediate_output.mfe'):
      if (("." in line) or ('(('in line) or ('))'in line)):        
        g=line.split()
        potential_structure=g[0]
        length_potential_structure=len(potential_structure)
        if length_potential_structure==sequence_length:
           structure=potential_structure    
    return(structure)

def determine_minimal_folding_energy(sequence, material):
    out_file=open('intermediate_output/intermediate_output.in', "w")
    out_file.write(standard_sequence_format(sequence))
    out_file.write("\n")
    out_file.close()
    
    st = os.stat('intermediate_output/intermediate_output.in')
    os.chmod('intermediate_output/intermediate_output.in', st.st_mode | stat.S_IEXEC)
    
    inp_file='intermediate_output/intermediate_output'
    bash_command='mfe -material '+material+' ' + inp_file
    subprocess.call(bash_command, shell=True)
    min_folding_energy=0
    for line in open('intermediate_output/intermediate_output.mfe'):    
        if (("%" not in line) and ("." in line)and ("-" in line)):
          h=line.split()
          folding_energy=float(h[0])   
          min_folding_energy=folding_energy
    
    return(min_folding_energy)

def determine_probability_secondary_structure(sequence, structure, material):
    out_file=open('intermediate_output/intermediate_output.in', "w")
    out_file.write(standard_sequence_format(sequence))
    out_file.write("\n")
    out_file.write(structure)
    out_file.write("\n")
    out_file.close()
    
    st = os.stat('intermediate_output/intermediate_output.in')
    os.chmod('intermediate_output/intermediate_output.in', st.st_mode | stat.S_IEXEC)
    
    inp_file='intermediate_output/intermediate_output'
    bash_command='prob -material '+material+' ' + inp_file +"> intermediate_output/intermediate_output.prob"
    subprocess.call(bash_command, shell=True)
    
    for line in open('intermediate_output/intermediate_output.prob'):
      if '%' not in line:
          f=line.split()
          str_probability=f[0]
          intermediart_probability=str_probability.replace(',', '.')
          probability=float(intermediart_probability)
          
    return(probability)
   
def get_complex_free_energy(iSBH_sequence, trigger_sequence, material):
    out_file=open('intermediate_output/intermediate_output_t.in', "w")
    out_file.write('2 ')
    out_file.write("\n")
    out_file.write(standard_sequence_format(iSBH_sequence))
    out_file.write("\n")
    out_file.write(standard_sequence_format(trigger_sequence))
    out_file.write("\n")
    out_file.write("1 2")
    out_file.write("\n")
    out_file.close()
    
    st = os.stat('intermediate_output/intermediate_output_t.in')
    os.chmod('intermediate_output/intermediate_output_t.in', st.st_mode | stat.S_IEXEC)
    
    inp_file='intermediate_output/intermediate_output_t'
    bash_command='mfe -multi -material '+material+' ' + inp_file 
    subprocess.call(bash_command, shell=True)
    
    for line in open('intermediate_output/intermediate_output_t.mfe'):
      if '-' in line:
          f=line.split()
          complex_energy=f[0].replace('\r', '')
    return(float(complex_energy))
    
def count_mismatches(observed_seq, expected_seq):
    length_observed_sequence=len(observed_seq)
    length_expected_sequence=len(expected_seq)
    if length_observed_sequence==length_expected_sequence:
        index=0
        mismatch_count=0
        list_mismatches=[0,0,0,0]
        while index < length_observed_sequence:
            if observed_seq[index]==expected_seq[index]:
                mismatch_count=mismatch_count+0
            else:
                mismatch_count=mismatch_count+1
                if mismatch_count <5:
                  list_mismatches[mismatch_count-1]=index+1
            index=index+1
    else:
        mismatch_count="error"
    
    return(mismatch_count, list_mismatches[0], list_mismatches[1],list_mismatches[2],list_mismatches[3])
    
def return_iSBH_fold_specified_design(design):
    fold="check fold name- non-standard format"
    if design=='SL14_iSBH':
        fold="((((((((((..(((..(((..............)))..)))..))))))))))"
    if design=='SL20_iSBH':
        fold="((((((((((..(((..(((....................)))..)))..))))))))))"
    if design=='SL25_iSBH':
        fold="((((((((((..(((..(((.........................)))..)))..))))))))))"
    if ("giraffe" in design) and ('SL14' in design):
        fold="((((((((((..(((..(((..(((..(((..............)))..)))..)))..)))..))))))))))"
    if ("giraffe" in design) and ('SL20' in design):
        fold="((((((((((..(((..(((..(((..(((....................)))..)))..)))..)))..))))))))))"
    if ("giraffe" in design) and ('SL25' in design):
        fold="((((((((((..(((..(((..(((..(((.........................)))..)))..)))..)))..))))))))))"
    return fold

def attempt_sequence_optimisation(RNA_sequence, desired_fold, beginning_flexible_seq, end_flexible_seq):
    if ((int(beginning_flexible_seq) <= int(end_flexible_seq)) and ((end_flexible_seq-1)<len(RNA_sequence)) and (len(RNA_sequence)==len(desired_fold))):
      fold_input_sequence=determine_mfe_structure(RNA_sequence, 'rna')
      index_beginning=beginning_flexible_seq-1
      index_end=end_flexible_seq-1
      unmodified_sequence_part_1=RNA_sequence[0:index_beginning]
      sequence_to_be_modified=RNA_sequence[index_beginning:index_end]
      unmodified_sequence_part_2=RNA_sequence[index_end:]
    
      desired_sub_structure_fold=desired_fold[index_beginning:index_end]
      actual_sub_structure_fold=fold_input_sequence[index_beginning:index_end]
      
      list_MM_indexes=[]
      number_mismatches=count_mismatches(desired_sub_structure_fold, actual_sub_structure_fold)[0]      
      MM_1_position=count_mismatches(desired_sub_structure_fold, actual_sub_structure_fold)[1]
      if MM_1_position !=0:
          list_MM_indexes.append(MM_1_position)
      MM_2_position=count_mismatches(desired_sub_structure_fold, actual_sub_structure_fold)[2]
      if MM_2_position !=0:
          list_MM_indexes.append(MM_2_position)
      MM_3_position=count_mismatches(desired_sub_structure_fold, actual_sub_structure_fold)[3]
      if MM_3_position !=0:
          list_MM_indexes.append(MM_3_position)
      MM_4_position=count_mismatches(desired_sub_structure_fold, actual_sub_structure_fold)[4]
      if MM_4_position !=0:
          list_MM_indexes.append(MM_4_position)
     
      if number_mismatches == 0:
          optimised_sequence=RNA_sequence
          
      else:         
          index=0
          list_optimised_sequences=[]
          number_sequences_to_be_screened=4**(len(list_MM_indexes))
          while index< number_sequences_to_be_screened:
              modified_subsequence=get_mismatched_RNA_sequence(sequence_to_be_modified,MM_1_position, MM_2_position, MM_3_position, MM_4_position)
              attempt_optimise_sequence=unmodified_sequence_part_1+modified_subsequence+unmodified_sequence_part_2
              structure_modification_attempt=determine_mfe_structure(attempt_optimise_sequence, 'rna')
              if structure_modification_attempt==desired_fold:
                  list_optimised_sequences.append(attempt_optimise_sequence)
              index=index+1
              
              if len(list_optimised_sequences)>0:
                  optimised_sequence=list_optimised_sequences[0]
              else:
                  optimised_sequence="XXXXX"
      return optimised_sequence           
    else:
      return "Error- Check if sequence and fold have save length! Check beginning and end of sequence!"

def reverse_complement(seq):
   complementary_sequence_list=[]
   index=0
   while index < len(seq):
            character=seq[index]
            if character=="A":
                complementary_base="T"
            if character=="T":
                complementary_base="A"
            if character=="C":
                complementary_base="G"
            if character=="G":
                complementary_base="C"
            complementary_sequence_list.append(complementary_base)
            index=index+1
            complementary_sequence=""
        
            for element in complementary_sequence_list:
              complementary_sequence=complementary_sequence+str(element)
            complementary_reverse=complementary_sequence[::-1] 
   return complementary_reverse

def ISBH_split(fw_isbh, rw_isbh):
    #initial split of the sequences
    if len(fw_isbh)%2 == 0:
      f1=fw_isbh[0:(len(fw_isbh)/2)]
      f2=fw_isbh[(len(fw_isbh)/2):len(fw_isbh)]
    else:
      f1=fw_isbh[0:((len(fw_isbh)-1)/2)]
      f2=fw_isbh[((len(fw_isbh)-1)/2):len(fw_isbh)] 
    if len(rw_isbh)%2 == 0:
      r1=rw_isbh[0:(len(rw_isbh)/2)]
      r2=rw_isbh[(len(rw_isbh)/2):len(rw_isbh)]
    else:
      r1=rw_isbh[0:((len(rw_isbh)-1)/2)]
      r2=rw_isbh[((len(rw_isbh)-1)/2):len(rw_isbh)]
      
    initial_overhang_1=f1[0:4]
    initial_overhang_2=r1[0:4]
    initial_overhang_3=f2[0:4]
    initial_overhang_4=r2[0:4]
    
    if (initial_overhang_1 != initial_overhang_3) or (initial_overhang_2 != initial_overhang_4):
       return (f1, r2, f2, r1)
    else:
       #second split of the sequences   
       if len(fw_isbh)%2 == 0:
         f1_b=fw_isbh[0:((len(fw_isbh)/2)-1)]
         f2_b=fw_isbh[((len(fw_isbh)/2)-1):len(fw_isbh)]
       else:
         f1_b=fw_isbh[0:(((len(fw_isbh)-1)/2)-1)]
         f2_b=fw_isbh[(((len(fw_isbh)-1)/2)-1):len(fw_isbh)]
       
       if len(rw_isbh)%2 == 0:
         r1_b=rw_isbh[0:((len(rw_isbh)/2)-1)]
         r2_b=rw_isbh[((len(rw_isbh)/2)-1):len(rw_isbh)]
       else:
         r1_b=rw_isbh[0:(((len(rw_isbh)-1)/2)-1)]
         r2_b=rw_isbh[(((len(rw_isbh)-1)/2)-1):len(rw_isbh)]
         
       second_overhang_1=f1_b[0:4]
       second_overhang_2=r1_b[0:4]
       second_overhang_3=f2_b[0:4]
       second_overhang_4=r2_b[0:4]
       if (second_overhang_1 != second_overhang_3) or (second_overhang_2 != second_overhang_4):        
         return (f1_b, r2_b, f2_b, r1_b)
       else: 
         #third split of the sequences           
         if len(fw_isbh)%2 == 0:
           f1_c=fw_isbh[0:((len(fw_isbh)/2)-2)]
           f2_c=fw_isbh[((len(fw_isbh)/2)-2):len(fw_isbh)]
         else:
           f1_c=fw_isbh[0:(((len(fw_isbh)-1)/2)-2)]
           f2_c=fw_isbh[(((len(fw_isbh)-1)/2)-2):len(fw_isbh)]
           
         if len(rw_isbh)%2 == 0:
           r1_c=rw_isbh[0:((len(rw_isbh)/2)-2)]
           r2_c=rw_isbh[((len(rw_isbh)/2)-2):len(rw_isbh)]
         else:
           r1_c=rw_isbh[0:(((len(rw_isbh)-1)/2)-2)]
           r2_c=rw_isbh[(((len(rw_isbh)-1)/2)-2):len(rw_isbh)]         
         third_overhang_1=f1_c[0:4]
         third_overhang_2=r1_c[0:4]
         third_overhang_3=f2_c[0:4]
         third_overhang_4=r2_c[0:4]
         if (third_overhang_1 != third_overhang_3) or (third_overhang_2 != third_overhang_4):        
           return (f1_c, r2_c, f2_c, r1_c)
         else:
           return ("error","error","error","error") 

def generate_plasmid_map(name_plasmid_backbone, restriction_enzyme_a, restriction_enzyme_b, FW_oligo, RW_oligo):
  for line in open('inputs/digested_plasmid_sequences.txt').readlines():
    key_search=name_plasmid_backbone+' '+restriction_enzyme_a+' '+restriction_enzyme_b
    if key_search in line:
      f=line.split()
      fragment_1_fw=standard_sequence_format(f[3])
      fragment_2_fw=standard_sequence_format(f[4])
      fragment_1_rw=standard_sequence_format(f[5])
      fragment_2_rw=standard_sequence_format(f[6])

      assembled_top_strand=fragment_1_fw+FW_oligo+fragment_2_fw
      assembled_bottom_strand=reverse_complement(fragment_1_rw+reverse_sequence(RW_oligo)+fragment_2_rw)
      
      
      if assembled_top_strand == reverse_sequence(assembled_bottom_strand):
          assembled_map=assembled_top_strand
          return(assembled_map)
      else:
          print('!!! Error generating plasmid map')

def find_DNA_sequence_seq_or_dna_files(filename):
    for line in open(filename).readlines():
      if ('A' in line) and ('C' in line) and ('T' in line) and ('G' in line) and ('EF' not in line):
        DNA_sequence=standard_sequence_format(line)
    return(DNA_sequence)

def analyse_seq_results(plasmid_map, sequencing_result, insert_seq):
    plasmid_sequence=find_DNA_sequence_seq_or_dna_files(plasmid_map)
    sequencing_sequence=find_DNA_sequence_seq_or_dna_files(sequencing_result)
    score=0
    if insert_seq in plasmid_sequence:
      score=score+1
      start_position_insert_in_plasmid_map=int(plasmid_sequence.index(insert_seq))
      index_before=start_position_insert_in_plasmid_map-30
      index_after_end=start_position_insert_in_plasmid_map+len(insert_seq)+500
      longer_sequence_to_search=plasmid_sequence[index_before: index_after_end]
      if longer_sequence_to_search in sequencing_sequence:
         score=score+1
      if len(sequencing_sequence) > 500:
        score=score+1
      if insert_seq in sequencing_sequence:
        score=score+1
      if score == 4:
        result="Sequencing test passed"
      elif score==3:
        result="!!! MANUALLY CHECK SEQUENCE"
      else:
        result="!!!!!!!!!! SEQUENCING TEST HAS FAILED"
    return(result)

def splitting_trigger_sequence(name, sequence):
    trigger=standard_sequence_format(sequence)
    if (("SL14" in name) and ("original_design" in name)):
        if len(sequence) == 34:
            loop=trigger[0:14]
            giraffe_extension=''
            spacer_star=trigger[14:]
        if len(sequence) == 33:
            loop=trigger[0:14]
            giraffe_extension=''
            spacer_star=trigger[14:]
    if (("SL20" in name) and ("original_design" in name)):
        if len(sequence) == 40:   
            loop=trigger[0:20]
            giraffe_extension=''
            spacer_star=trigger[20:]

    if (("SL25" in name) and ("original_design" in name)):
        if len(sequence) == 45: 
            loop=trigger[0:25]
            giraffe_extension=''
            spacer_star=trigger[25:]

    if ("SL14" in name) and ("giraffe" in name):
        if len(sequence) == 44:
            loop=trigger[0:14]
            giraffe_extension=trigger[14:24]
            spacer_star=trigger[24:]
         
    if ("SL20" in name) and ("giraffe" in name):
        if len(sequence) == 50:   
            loop=trigger[0:20]
            giraffe_extension=trigger[20:30]
            spacer_star=trigger[30:]
   
    if ("SL25" in name) and ("giraffe" in name):
        if len(sequence) == 55: 
            loop=trigger[0:25]
            giraffe_extension=trigger[25:35]
            spacer_star=trigger[35:]
    if (loop+giraffe_extension+spacer_star)==trigger:
        return(loop, giraffe_extension, spacer_star)
    else:
        return ('error_check_splitt_trigger_function','error','error')

def random_nucleotide_sequence_except_input(initial_sequence, dna_or_rna):
   if dna_or_rna=='dna':
     index=0
     random_sequence=''
     while index < len(initial_sequence):
       initial_nucleotide=initial_sequence[index]
       if initial_nucleotide=='A':
          random_nucleotide=random.choice(['C', 'T', 'G'])
       if initial_nucleotide=='C':
          random_nucleotide=random.choice(['A', 'T', 'G'])
       if initial_nucleotide=='T':
          random_nucleotide=random.choice(['C', 'A', 'G'])
       if initial_nucleotide=='G':
          random_nucleotide=random.choice(['C', 'T', 'A']) 
       random_sequence=random_sequence+random_nucleotide
       index=index+1
   if dna_or_rna=='rna':
     index=0
     random_sequence=''
     while index < len(initial_sequence):
       initial_nucleotide=initial_sequence[index]
       if initial_nucleotide=='A':
          random_nucleotide=random.choice(['C', 'U', 'G'])
       if initial_nucleotide=='C':
          random_nucleotide=random.choice(['A', 'U', 'G'])
       if initial_nucleotide=='U':
          random_nucleotide=random.choice(['C', 'A', 'G'])
       if initial_nucleotide=='G':
          random_nucleotide=random.choice(['C', 'U', 'A'])  
       random_sequence=random_sequence+random_nucleotide
       index=index+1
   return (random_sequence)

def return_random_sequence_exccept_desired_matches(any_sequence, dna_or_rna, M1, M2, M3, M4):  #always need to add 4 matches, if need less than 4 matches, repeat a number twice
    sequence=standard_sequence_format(any_sequence)
    list_of_matched_indexes=[]
    if M1 !=0:
        M1_index=M1-1
        list_of_matched_indexes.append(M1_index)
    if M2 !=0:
        M2_index=M2-1
        list_of_matched_indexes.append(M2_index)
    if M3 !=0:
        M3_index=M3-1
        list_of_matched_indexes.append(M3_index)
    if M3 !=0:
        M4_index=M4-1
        list_of_matched_indexes.append(M4_index)
    index=0
    desired_sequence=''
    while index < len(sequence):
        initial_nucleotide=sequence[index]
        if (index==list_of_matched_indexes[0]) or (index==list_of_matched_indexes[1]) or (index==list_of_matched_indexes[2]) or (index==list_of_matched_indexes[3]):
            desired_nucleotide=initial_nucleotide
        else:
            if dna_or_rna =='dna':
                desired_nucleotide=random_nucleotide_sequence_except_input(initial_nucleotide, 'dna').lower()
            if dna_or_rna =='rna':
                desired_nucleotide=random_nucleotide_sequence_except_input(initial_nucleotide, 'rna').lower()
        desired_sequence=desired_sequence+desired_nucleotide
        index=index+1
    return(desired_sequence)

def add_BbsI_restriction_sites_and_extra_G_U6_sg(FW_DNA_seq): #for iSBHs, an extra C complementary to tracrRNA is not added
    if FW_DNA_seq[0]=='G':
        U6_G=''
        U6_complementary_C=''
    else:
        U6_G='G'
        U6_complementary_C='C'
    FW_sequence_with_restriction_sites="CACC"+U6_G+FW_DNA_seq
    RW_sequence_with_restriction_sites="AAAC"+get_DNA_reverse_complement(FW_DNA_seq)+U6_complementary_C
    return(FW_sequence_with_restriction_sites, RW_sequence_with_restriction_sites)

def add_BbsI_restriction_sites_TM2(FW_DNA_seq): #an extra U6 G is not added
    FW_sequence_with_restriction_sites="CCAC"+FW_DNA_seq
    RW_sequence_with_restriction_sites="GCAT"+get_DNA_reverse_complement(FW_DNA_seq)
    return(FW_sequence_with_restriction_sites, RW_sequence_with_restriction_sites)

def add_XbaI_and_AscI_restriction_sites_1xCTS(FW_CTS_sequence):
    if len(FW_CTS_sequence)!=23:
        print '!!! Potential error- no conventional 23 nt CTS'
    FW_sequence_with_restriction_sites="CTAGA"+FW_CTS_sequence+"GG"
    RW_sequence_with_restriction_sites="CGCGCC"+get_DNA_reverse_complement(FW_CTS_sequence)+"T"
    return(FW_sequence_with_restriction_sites, RW_sequence_with_restriction_sites)

def return_number_matches_2_sequences(sequence1, sequence2):
  if len(sequence1)== len(sequence2):
      sequence_1=standard_sequence_format(sequence1)
      sequence_2=standard_sequence_format(sequence2)
      index=0
      number_matches=0
      while index<len(sequence_1):
          if sequence_1[index]==sequence_2[index]:
              number_matches=number_matches+1
          index=index+1
      return(number_matches)          
  else:
      return('error')

def return_percentage_similarity_2_sequences(sequence1, sequence2):
      if len(sequence1)== len(sequence2):
        sequence_length=float(len(sequence1))
        number_matches=float(return_number_matches_2_sequences(sequence1, sequence2))
        percentage_similarity=(number_matches/sequence_length)*100.00
        return(percentage_similarity)
      else:
         return ('error')

def calculate_Tm(primer):
  primer_sequence=standard_sequence_format(primer)
  A_count=primer_sequence.count('A')
  T_count=primer_sequence.count('C')
  C_count=primer_sequence.count('T')
  G_count=primer_sequence.count('G')
  Tm=float((2*(A_count+T_count))+(4*(C_count+G_count)))
  return(Tm)     
  
def find_nearest(array, value):
    array=np.asarray(array)
    idx=(np.abs(array-value)).argmin()
    return (array[idx])

def design_annealing_primers(given_sequence):
    FW_sequence=standard_sequence_format(given_sequence)
    RW_sequence=get_DNA_reverse_complement(FW_sequence)
    min_Tm=50.00
    max_Tm=72.00
    if len(given_sequence) >=50:
       index=17
       FW_primer_seq_list=[]
       FW_primer_Tms=[]
       while index <= 30:
           if (FW_sequence[index]=="G") or (FW_sequence[index]=="C"):
               primer=FW_sequence[0:(index+1)]
               Tm_content=float(calculate_Tm(primer))
               if (min_Tm <= Tm_content) and (Tm_content <=max_Tm):
                   FW_primer_seq_list.append(primer)
                   FW_primer_Tms.append(Tm_content)
           index=index+1
       
       index_2=17
       RW_primer_seq_list=[]
       RW_primer_Tms=[]
       while index_2 <= 30:
           if (RW_sequence[index_2]=="G") or (RW_sequence[index_2]=="C"):
               primer_2=RW_sequence[0:(index_2+1)]
               Tm_content=float(calculate_Tm(primer_2))
               if (min_Tm <= Tm_content) and (Tm_content <=max_Tm):
                   RW_primer_seq_list.append(primer_2)
                   RW_primer_Tms.append(Tm_content)
           index_2=index_2+1  
           
       Tm_FW=np.array(FW_primer_Tms)
       Tm_RW=np.array(RW_primer_Tms)
       average_Tm_FW_primers=np.average(Tm_FW)
       average_Tm_RW_primers=np.average(Tm_RW)
       average_Tms_all_primers=(average_Tm_FW_primers+average_Tm_RW_primers)/2
       
       FW_primer_optimal_Tm=find_nearest(Tm_FW, average_Tms_all_primers)
       RW_primer_optimal_Tm=find_nearest(Tm_RW, average_Tms_all_primers)
       
       id_FW_prim=FW_primer_Tms.index(FW_primer_optimal_Tm)
       id_RW_prim=RW_primer_Tms.index(RW_primer_optimal_Tm)
       
       FW_primer=FW_primer_seq_list[id_FW_prim]
       RW_primer=RW_primer_seq_list[id_RW_prim]
       
    else:
        return ('error')
    return(FW_primer, RW_primer)

def return_list_non_identical_nucleotides(DNA_nucleotide):
    nucleotide=standard_sequence_format(DNA_nucleotide)
    if nucleotide=='A':
        nucleotide_list=['C','T', 'G']
    if nucleotide=='C':
        nucleotide_list=['A','T', 'G']
    if nucleotide=='T':
        nucleotide_list=['A','C', 'G']
    if nucleotide=='G':
        nucleotide_list=['A','C','T']
    return(nucleotide_list)
  

def return_all_bulge_combinations_except_input(two_nt_bulge):
    if len(two_nt_bulge)==2:
        first_nucleotide=two_nt_bulge[0]
        second_nucleotide=two_nt_bulge[1]
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                new_bulge=nucleotide_1+nucleotide_2
                list_all_potential_combinations.append(new_bulge)
        return(list_all_potential_combinations)
    else:
        return('error')  
        
def return_potential_giraffe_sequence(potential_trigger_full_sequence, desired_scaffold):
    if 'P1' in desired_scaffold:
        spacer_star='TGCTTCGCTAGTCGC'
        spacer='TCGCGTGTAGCGAAGCA'    
    if 'P2' in desired_scaffold:
        spacer_star='AGGACAGTACAGCGA'
        spacer='AGTCGGAGTACTGTCCT'
    if 'G2' in desired_scaffold:
        spacer_star='TCCTGCTAGTACTTA'
        spacer='GGTAATGACTAGCAGGA'  
    if 'G3' in desired_scaffold:
        spacer_star='GAGGCGGGAAGTATG'
        spacer='GACATTATTCCCGCCTC'
    if 'G5' in desired_scaffold:
        spacer_star='GACTCTTACTCCAGT'
        spacer='CTACTTCAGTAAGAGTC' 
    
    SL14_fold_giraffe_fold=return_iSBH_fold_specified_design('SL14_giraffe_design')
    
    #building giraffe designs without bulges
    trigger_principal_sensing_component=potential_trigger_full_sequence[0:29]
    sensing_module=get_DNA_reverse_complement(trigger_principal_sensing_component)
    giraffe_extension_star=sensing_module[5:15]
    three_nt_spacer_component=get_DNA_reverse_complement(sensing_module[2:5])   
    giraffe_extension=get_DNA_reverse_complement(giraffe_extension_star)
    
    giraff_iSBH_no_bulges=spacer_star+sensing_module+giraffe_extension+three_nt_spacer_component+spacer
    
    #systematic generation of all possible bulge nucleotide combinations
    nucleotides_first_bulge=giraff_iSBH_no_bulges[47:49]
    list_all_potential_combinations_first_bulge=return_all_bulge_combinations_except_input(nucleotides_first_bulge)
    
    nucleotides_second_bulge=giraff_iSBH_no_bulges[52:54]
    list_all_potential_combinations_second_bulge=return_all_bulge_combinations_except_input(nucleotides_second_bulge)
    
    #generating a list of giraffe folds with all possible bulge combinations
    list_all_potential_designs_no_extra_trigger_complementarity=[]
    for bulge_1 in list_all_potential_combinations_first_bulge:
        for bulge_2 in list_all_potential_combinations_second_bulge:
            potential_iSBH_giraffe_design=giraff_iSBH_no_bulges[0:47]+bulge_1+giraff_iSBH_no_bulges[49:52]+bulge_2+giraff_iSBH_no_bulges[54:]
            list_all_potential_designs_no_extra_trigger_complementarity.append(potential_iSBH_giraffe_design)
    
    #for list all possible bulge combinations, also generate iSBHs with the second spacer * bulge complementary to the trigger
    list_all_potential_designs_extra_trigger_complementarity=[]
    trigger_complementary_second_spacer_bulge=get_DNA_reverse_complement(potential_trigger_full_sequence[32:34])
    
    for potential_iSBH in list_all_potential_designs_no_extra_trigger_complementarity:
        extra_complementarity_iSBH=potential_iSBH[0:10]+trigger_complementary_second_spacer_bulge+potential_iSBH[12:]
        list_all_potential_designs_extra_trigger_complementarity.append(extra_complementarity_iSBH)

    #generate a list comprising all iSBH designs that fold into the right configuration
    list_all_designs_generated=list_all_potential_designs_no_extra_trigger_complementarity+list_all_potential_designs_extra_trigger_complementarity
    list_all_designs_right_fold=[]
    for every_sequence in list_all_designs_generated:
        sequence_structure=determine_mfe_structure(every_sequence, 'rna')
        if sequence_structure==SL14_fold_giraffe_fold:
            list_all_designs_right_fold.append(every_sequence)
    
    return(list_all_designs_right_fold)

def return_potential_giraffe_sequence_endogenous_spacers(potential_trigger_44nt_full_sequence, desired_17_nt_spacer, desired_15_nt_spacer_star_with_bulges):
    spacer_star=desired_15_nt_spacer_star_with_bulges
    spacer=desired_17_nt_spacer 
    
    SL14_fold_giraffe_fold=return_iSBH_fold_specified_design('SL14_giraffe_design')
    
    #building giraffe designs without bulges
    trigger_principal_sensing_component=potential_trigger_44nt_full_sequence[0:29]
    sensing_module=get_DNA_reverse_complement(trigger_principal_sensing_component)
    giraffe_extension_star=sensing_module[5:15]
    three_nt_spacer_component=get_DNA_reverse_complement(sensing_module[2:5])   
    giraffe_extension=get_DNA_reverse_complement(giraffe_extension_star)
    
    giraff_iSBH_no_bulges=spacer_star+sensing_module+giraffe_extension+three_nt_spacer_component+spacer
    
    #systematic generation of all possible bulge nucleotide combinations
    nucleotides_first_bulge=giraff_iSBH_no_bulges[47:49]
    list_all_potential_combinations_first_bulge=return_all_bulge_combinations_except_input(nucleotides_first_bulge)
    
    nucleotides_second_bulge=giraff_iSBH_no_bulges[52:54]
    list_all_potential_combinations_second_bulge=return_all_bulge_combinations_except_input(nucleotides_second_bulge)
    
    #generating a list of giraffe folds with all possible bulge combinations
    list_all_potential_designs_no_extra_trigger_complementarity=[]
    for bulge_1 in list_all_potential_combinations_first_bulge:
        for bulge_2 in list_all_potential_combinations_second_bulge:
            potential_iSBH_giraffe_design=giraff_iSBH_no_bulges[0:47]+bulge_1+giraff_iSBH_no_bulges[49:52]+bulge_2+giraff_iSBH_no_bulges[54:]
            list_all_potential_designs_no_extra_trigger_complementarity.append(potential_iSBH_giraffe_design)
    
    #for list all possible bulge combinations, also generate iSBHs with the second spacer * bulge complementary to the trigger
    list_all_potential_designs_extra_trigger_complementarity=[]
    trigger_complementary_second_spacer_bulge=get_DNA_reverse_complement(potential_trigger_44nt_full_sequence[32:34])
    
    for potential_iSBH in list_all_potential_designs_no_extra_trigger_complementarity:
        extra_complementarity_iSBH=potential_iSBH[0:10]+trigger_complementary_second_spacer_bulge+potential_iSBH[12:]
        list_all_potential_designs_extra_trigger_complementarity.append(extra_complementarity_iSBH)

    #generate a list comprising all iSBH designs that fold into the right configuration
    list_all_designs_generated=list_all_potential_designs_no_extra_trigger_complementarity+list_all_potential_designs_extra_trigger_complementarity
    list_all_designs_right_fold=[]
    for every_sequence in list_all_designs_generated:
        sequence_structure=determine_mfe_structure(every_sequence, 'rna')
        if sequence_structure==SL14_fold_giraffe_fold:
            list_all_designs_right_fold.append(every_sequence)
    
    return(list_all_designs_right_fold)

def return_potential_SL25_giraffe_sequence_endogenous_spacers(potential_trigger_55nt_full_sequence, desired_17_nt_spacer, desired_15_nt_spacer_star_with_bulges):
    spacer_star=desired_15_nt_spacer_star_with_bulges
    spacer=desired_17_nt_spacer 
    
    SL14_fold_giraffe_fold=return_iSBH_fold_specified_design('SL25_giraffe_design')
    
    #building giraffe designs without bulges
    trigger_principal_sensing_component=potential_trigger_55nt_full_sequence[0:40]
    sensing_module=get_DNA_reverse_complement(trigger_principal_sensing_component)
    giraffe_extension_star=sensing_module[5:15]
    three_nt_spacer_component=get_DNA_reverse_complement(sensing_module[2:5])   
    giraffe_extension=get_DNA_reverse_complement(giraffe_extension_star)
    
    giraff_iSBH_no_bulges=spacer_star+sensing_module+giraffe_extension+three_nt_spacer_component+spacer
    
    #systematic generation of all possible bulge nucleotide combinations
    nucleotides_first_bulge=giraff_iSBH_no_bulges[(47+11):(49+11)]
    list_all_potential_combinations_first_bulge=return_all_bulge_combinations_except_input(nucleotides_first_bulge)
    
    nucleotides_second_bulge=giraff_iSBH_no_bulges[(52+11):(54+11)]
    list_all_potential_combinations_second_bulge=return_all_bulge_combinations_except_input(nucleotides_second_bulge)
    
    #generating a list of giraffe folds with all possible bulge combinations
    list_all_potential_designs_no_extra_trigger_complementarity=[]
    for bulge_1 in list_all_potential_combinations_first_bulge:
        for bulge_2 in list_all_potential_combinations_second_bulge:
            potential_iSBH_giraffe_design=giraff_iSBH_no_bulges[0:(47+11)]+bulge_1+giraff_iSBH_no_bulges[(49+11):(52+11)]+bulge_2+giraff_iSBH_no_bulges[(54+11):]
            list_all_potential_designs_no_extra_trigger_complementarity.append(potential_iSBH_giraffe_design)
    
    #for list all possible bulge combinations, also generate iSBHs with the second spacer * bulge complementary to the trigger
    list_all_potential_designs_extra_trigger_complementarity=[]
    trigger_complementary_second_spacer_bulge=get_DNA_reverse_complement(potential_trigger_55nt_full_sequence[(32):(34)])
    
    for potential_iSBH in list_all_potential_designs_no_extra_trigger_complementarity:
        extra_complementarity_iSBH=potential_iSBH[0:10]+trigger_complementary_second_spacer_bulge+potential_iSBH[12:]
        list_all_potential_designs_extra_trigger_complementarity.append(extra_complementarity_iSBH)

    #generate a list comprising all iSBH designs that fold into the right configuration
    list_all_designs_generated=list_all_potential_designs_no_extra_trigger_complementarity+list_all_potential_designs_extra_trigger_complementarity
    list_all_designs_right_fold=[]
    for every_sequence in list_all_designs_generated:
        sequence_structure=determine_mfe_structure(every_sequence, 'rna')
        if sequence_structure==SL14_fold_giraffe_fold:
            list_all_designs_right_fold.append(every_sequence)
    
    return(list_all_designs_right_fold)

def return_potential_SL20_giraffe_sequence_endogenous_spacers(potential_trigger_50nt_full_sequence, desired_17_nt_spacer, desired_15_nt_spacer_star_with_bulges):
    spacer_star=desired_15_nt_spacer_star_with_bulges
    spacer=desired_17_nt_spacer 
    
    SL14_fold_giraffe_fold=return_iSBH_fold_specified_design('SL20_giraffe_design')
    
    #building giraffe designs without bulges
    trigger_principal_sensing_component=potential_trigger_50nt_full_sequence[0:35]
    sensing_module=get_DNA_reverse_complement(trigger_principal_sensing_component)
    giraffe_extension_star=sensing_module[5:15]
    three_nt_spacer_component=get_DNA_reverse_complement(sensing_module[2:5])   
    giraffe_extension=get_DNA_reverse_complement(giraffe_extension_star)
    
    giraff_iSBH_no_bulges=spacer_star+sensing_module+giraffe_extension+three_nt_spacer_component+spacer
    
    #systematic generation of all possible bulge nucleotide combinations
    nucleotides_first_bulge=giraff_iSBH_no_bulges[(47+6):(49+6)]
    list_all_potential_combinations_first_bulge=return_all_bulge_combinations_except_input(nucleotides_first_bulge)
    
    nucleotides_second_bulge=giraff_iSBH_no_bulges[(52+6):(54+6)]
    list_all_potential_combinations_second_bulge=return_all_bulge_combinations_except_input(nucleotides_second_bulge)
    
    #generating a list of giraffe folds with all possible bulge combinations
    list_all_potential_designs_no_extra_trigger_complementarity=[]
    for bulge_1 in list_all_potential_combinations_first_bulge:
        for bulge_2 in list_all_potential_combinations_second_bulge:
            potential_iSBH_giraffe_design=giraff_iSBH_no_bulges[0:(47+6)]+bulge_1+giraff_iSBH_no_bulges[(49+6):(52+6)]+bulge_2+giraff_iSBH_no_bulges[(54+6):]
            list_all_potential_designs_no_extra_trigger_complementarity.append(potential_iSBH_giraffe_design)
    
    #for list all possible bulge combinations, also generate iSBHs with the second spacer * bulge complementary to the trigger
    list_all_potential_designs_extra_trigger_complementarity=[]
    trigger_complementary_second_spacer_bulge=get_DNA_reverse_complement(potential_trigger_50nt_full_sequence[(32):(34)])
    
    for potential_iSBH in list_all_potential_designs_no_extra_trigger_complementarity:
        extra_complementarity_iSBH=potential_iSBH[0:10]+trigger_complementary_second_spacer_bulge+potential_iSBH[12:]
        list_all_potential_designs_extra_trigger_complementarity.append(extra_complementarity_iSBH)

    #generate a list comprising all iSBH designs that fold into the right configuration
    list_all_designs_generated=list_all_potential_designs_no_extra_trigger_complementarity+list_all_potential_designs_extra_trigger_complementarity
    list_all_designs_right_fold=[]
    for every_sequence in list_all_designs_generated:
        sequence_structure=determine_mfe_structure(every_sequence, 'rna')
        if sequence_structure==SL14_fold_giraffe_fold:
            list_all_designs_right_fold.append(every_sequence)
    
    return(list_all_designs_right_fold)

def return_all_sequence_combinations_desired_seq_up_to_4nt(seq_length):
    list_all_1_nt_combination=['A','C','T','G']
    list_all_2_nt_combination=[]
    list_all_3_nt_combination=[]
    list_all_4_nt_combination=[]
    
    for nt_1 in list_all_1_nt_combination:
        for nt_2 in list_all_1_nt_combination:
            list_all_2_nt_combination.append(nt_1+nt_2)
            for nt_3 in list_all_1_nt_combination:
                list_all_3_nt_combination.append(nt_1+nt_2+nt_3)
                for nt_4 in list_all_1_nt_combination:
                  list_all_4_nt_combination.append(nt_1+nt_2+nt_3+nt_4)
    if   seq_length==1:
       list_all_sequence_combinations=  list_all_1_nt_combination
       
    if   seq_length==2:
       list_all_sequence_combinations=  list_all_2_nt_combination 

    if   seq_length==3:
       list_all_sequence_combinations=  list_all_3_nt_combination
       
    if   seq_length==4:
       list_all_sequence_combinations=  list_all_4_nt_combination
       
    return (list_all_sequence_combinations)

######new functions
######new functions
######new functions

def generate_symbolic_iSBH_fold_from_given_iSBH_components_no_toehold(giraffe_ext_star_fold, giraffe_extension_fold,loop_length):
    
    spacer_star_fold='((((((((((..(((..((('
    spacer_star=')))..)))..))))))))))'
    
    loop_fold='.'*loop_length
    
    desired_fold=spacer_star_fold+giraffe_ext_star_fold+loop_fold+giraffe_extension_fold+spacer_star
   
    trigger_sub_sequence_length_excluding_toehold=len(giraffe_ext_star_fold)+loop_length
    
    length_iSBH=len(desired_fold)
    print desired_fold
    return (desired_fold, trigger_sub_sequence_length_excluding_toehold, length_iSBH, loop_length)

def generate_iSBH_symbolic_fold_from_given_inputs_with_toehold(length_base_pairing_tracr, giraffe_ext_star_full_fold, giraffe_extension_full_fold ,loop_length, toehold_length):
    
    extraG_fold='.'
    toehold_fold='.'*toehold_length    
    extra_interanctions_tracr='('*length_base_pairing_tracr
    
    spacer_star_fold='((((((((((..(((..((('
    spacer_star=')))..)))..))))))))))'
    
    loop_fold='.'*loop_length

    tracer_extra_interactions=')'*length_base_pairing_tracr
    
    desired_fold=extraG_fold+toehold_fold+extra_interanctions_tracr+spacer_star_fold+giraffe_ext_star_full_fold+loop_fold+giraffe_extension_full_fold+spacer_star+tracer_extra_interactions
   
    trigger_sub_sequence_length_excluding_toehold=len(giraffe_ext_star_full_fold)+loop_length
    
    length_iSBH=len(desired_fold)
    
    return (desired_fold, trigger_sub_sequence_length_excluding_toehold, length_iSBH, loop_length)

def generate_iSBH_sensing_modules_for_given_trigger_and_specifications(trigger_sequence, trigger_sub_sequence_length_excluding_toehold):
    list_sensing_modules=[]
    index=0
    while index <(len(trigger_sequence)-trigger_sub_sequence_length_excluding_toehold):
        potential_sensing_module=get_DNA_reverse_complement(trigger_sequence[index:(index+trigger_sub_sequence_length_excluding_toehold)])
        list_sensing_modules.append(potential_sensing_module)
        index=index+1
        
    return (list_sensing_modules)


def determine_most_stable_backfold_for_given_spacer(spacer_DNA_sequence):
    
    #determine minimal iSBH seq with a small loop and no bulges
    reverse_complement_spacer=get_DNA_reverse_complement(spacer_DNA_sequence)
    minimal_loop='GAAA'#'GAAA'
    minimal_iSBH_no_bulges=reverse_complement_spacer+minimal_loop+spacer_DNA_sequence

    #determine all potential nucleotide mis-matches that might lead to the desired fold
    
    desired_fold_minimal_iSBH=generate_symbolic_iSBH_fold_from_given_iSBH_components_no_toehold('', '',4)[0]
    first_two_nt_to_mis_match=minimal_iSBH_no_bulges[10:12]
    last_two_nt_to_mis_match=minimal_iSBH_no_bulges[15:17]
    
    list_all_potential_first_bulges=return_all_bulge_combinations_except_input(first_two_nt_to_mis_match)
    list_all_potential_second_bulges=return_all_bulge_combinations_except_input(last_two_nt_to_mis_match)
            
            #determine all potential bulge combinations and output assembled minimal iSBH
    
    list_all_potential_combinations_minimal_iSBHs_with_bulges=[]
    list_folding_probabilities_for_minimal_iSBHs=[]
    
    for bulge_1 in list_all_potential_first_bulges:
        first_two_nt_to_mis_match=bulge_1
        for bulge_2 in list_all_potential_second_bulges:
            last_two_nt_to_mis_match=bulge_2
            assembled_iSBH= minimal_iSBH_no_bulges[:10]+first_two_nt_to_mis_match+minimal_iSBH_no_bulges[12:15]+last_two_nt_to_mis_match+minimal_iSBH_no_bulges[17:]
            fold_assembled_iSBH=determine_mfe_structure(assembled_iSBH, 'rna')
            if fold_assembled_iSBH== desired_fold_minimal_iSBH:
                list_all_potential_combinations_minimal_iSBHs_with_bulges.append(assembled_iSBH)
                folding_probability=determine_probability_secondary_structure(assembled_iSBH, desired_fold_minimal_iSBH, 'rna')
                list_folding_probabilities_for_minimal_iSBHs.append(folding_probability)
   
    #pick best sequence
    index_best_fold=list_folding_probabilities_for_minimal_iSBHs.index(max(list_folding_probabilities_for_minimal_iSBHs))
    small_iSBH_sequence=list_all_potential_combinations_minimal_iSBHs_with_bulges[index_best_fold]
    spacer_star_DNA_sequence=small_iSBH_sequence[:20]
    return (spacer_star_DNA_sequence, spacer_DNA_sequence, small_iSBH_sequence)
    
def assemble_preliminary_iSBH_sequence_no_bulges(list_sensing_modules, spacer_star_DNA_sequence, spacer_DNA_sequence, desired_fold, folding_probability_cutting_threshold):
  spacer_star_fold='((((((((((..(((..((('
  spacer_star=')))..)))..))))))))))'
  loop_start=0
  loop_end=0
  index=0
  while index < len(desired_fold):
      if (desired_fold[(index-1)]=='(') and (desired_fold[index:(index+14)]=='..............'):
        loop_start=index
      if (desired_fold[index:(index+14)]=='..............') and (desired_fold[(index+15)]==')'):
        loop_end=index+14
      index=index+1
  
  loop_length=len(desired_fold[loop_start : loop_end])
  stem_length=len(desired_fold[ : loop_start])
  desired_provisional_fold=spacer_star_fold+((stem_length-20)*"(")+("."*loop_length)+(")"*(stem_length-20))+spacer_star
  list_potential_iSBHs_no_bulges=[]
    
  for sensing_module in list_sensing_modules:
       
        stem_star_component=sensing_module[0:-loop_length]
        stem_component=get_DNA_reverse_complement(stem_star_component)
        loop_component=sensing_module[-loop_length:]

        assembled_sequence=spacer_star_DNA_sequence+stem_star_component+loop_component+stem_component+spacer_DNA_sequence
        
        folding_assembled_sequence=determine_mfe_structure(assembled_sequence, 'rna')
        if desired_provisional_fold == folding_assembled_sequence:
          folding_probability=determine_probability_secondary_structure(assembled_sequence,desired_provisional_fold,'rna')
          if folding_probability > folding_probability_cutting_threshold:
              list_potential_iSBHs_no_bulges.append(assembled_sequence)
  return(list_potential_iSBHs_no_bulges)

def identifying_corresponding_MM_positions(desired_fold):
  list_MMs=[]
  loop_start=0
  loop_end=0
  index=0
  while index < len(desired_fold):
      if (desired_fold[(index-1)]=='(') and (desired_fold[index:(index+14)]=='..............'):
        loop_start=index
      if (desired_fold[index:(index+14)]=='..............') and (desired_fold[(index+15)]==')'):
        loop_end=index+14
      index=index+1
  
  stem_length=len(desired_fold[ : loop_start])
 
  list_MM_giraffe_star=[]
  index_2=18
  while index_2 < stem_length:

      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+2)]=='.(':
          MM=[index_2]
          list_MM_giraffe_star.append(MM)
          
      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+3)]=='..(':
          MM=[index_2,index_2+1]
          list_MM_giraffe_star.append(MM)   
          
      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+4)]=='...(':
          MM=[index_2,index_2+1,index_2+2]
          list_MM_giraffe_star.append(MM) 
          
      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+5)]=='....(':
          MM=[index_2,index_2+1,index_2+2,index_2+3]
          list_MM_giraffe_star.append(MM)          

      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+6)]=='.....(':
          MM=[index_2,index_2+1,index_2+2,index_2+3, index_2+4]
          list_MM_giraffe_star.append(MM) 

      if (desired_fold[index_2-1]=='(') and desired_fold[index_2:(index_2+7)]=='......(':
          MM=[index_2,index_2+1,index_2+2,index_2+3, index_2+4,index_2+5]
          list_MM_giraffe_star.append(MM) 
      
      index_2=index_2+1
  
  list_MM_giraffe=[]
  list_inverted_MM_giraffe=[]
  
  index_3=0  
  while index_3 < len(desired_fold)-18:
      if (desired_fold[index_3-1]==')') and desired_fold[index_3:(index_3+2)]=='.)':
          MMg=[index_3]
          list_MM_giraffe.append(MMg)
      
      if (desired_fold[index_3-1]==')') and desired_fold[index_3:(index_3+3)]=='..)':
          MMg=[index_3, index_3+1]
          list_MM_giraffe.append(MMg)
          
      if (desired_fold[index_3-1]==')') and desired_fold[index_3:(index_3+4)]=='...)':
          MMg=[index_3, index_3+1, index_3+2]
          list_MM_giraffe.append(MMg)       

      if (desired_fold[index_3-1]==')') and desired_fold[index_3:(index_3+5)]=='....)':
          MMg=[index_3, index_3+1, index_3+2, index_3+3]
          list_MM_giraffe.append(MMg)   
          
      if (desired_fold[index_3-1]==')') and desired_fold[index_3:(index_3+6)]=='.....)':
          MMg=[index_3, index_3+1, index_3+2, index_3+3, index_3+4]
          list_MM_giraffe.append(MMg)   
                    
          
      index_3=index_3+1  
      
  list_inverted_MM_giraffe=[] 
  index_4=len(list_MM_giraffe)-1
  while index_4 >= 0:
        list_inverted_MM_giraffe.append(list_MM_giraffe[index_4])
        index_4=index_4-1
  
  number_asymetrical_MMs=0
  
  index_5=len(list_MM_giraffe)-1
  while index_5>=0:
      list_corresponding_MMs=[list_MM_giraffe_star[index_5], list_inverted_MM_giraffe[index_5]]
      no_MM_giraffe_star=len(list_MM_giraffe_star[index_5])
      no_MM_giraffe=len(list_inverted_MM_giraffe[index_5])
      number_asymmetrical_bp_for_this_MM=no_MM_giraffe_star-no_MM_giraffe
      number_asymetrical_MMs=number_asymetrical_MMs+number_asymmetrical_bp_for_this_MM
      list_MMs.append(list_corresponding_MMs)
      index_5=index_5-1
      
  return(list_MMs, number_asymetrical_MMs)

def return_all_bulge_combinations_except_input_for_up_to_6nt_seq(bulge):
    
    #bulge=1nt
    if len(bulge)==1:
        first_nucleotide=bulge[0]
        list_all_potential_combinations=return_list_non_identical_nucleotides(first_nucleotide)
    
    #bulge=2nt
    if len(bulge)==2:
        first_nucleotide=bulge[0]
        second_nucleotide=bulge[1]
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                new_bulge=nucleotide_1+nucleotide_2
                list_all_potential_combinations.append(new_bulge)
    
    #bulge=3nt           
    if len(bulge)==3:
        first_nucleotide=bulge[0]
        second_nucleotide=bulge[1]
        third_nucleotide=bulge[2]
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_potential_nt_third_nt=return_list_non_identical_nucleotides(third_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                for nucleotide_3 in list_potential_nt_third_nt:
                  new_bulge=nucleotide_1+nucleotide_2+nucleotide_3
                  list_all_potential_combinations.append(new_bulge)          
    #bulge=4nt           
    if len(bulge)==4:
        first_nucleotide=bulge[0]
        second_nucleotide=bulge[1]
        third_nucleotide=bulge[2]
        fourth_nucleotide=bulge[3]
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_potential_nt_third_nt=return_list_non_identical_nucleotides(third_nucleotide)
        list_potential_nt_fourth_nt=return_list_non_identical_nucleotides(fourth_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                for nucleotide_3 in list_potential_nt_third_nt:
                    for nucleotide_4 in list_potential_nt_fourth_nt:
                       new_bulge=nucleotide_1+nucleotide_2+nucleotide_3+nucleotide_4
                       list_all_potential_combinations.append(new_bulge)                   
    #bulge=5nt           
    if len(bulge)==5:
        first_nucleotide=bulge[0]
        second_nucleotide=bulge[1]
        third_nucleotide=bulge[2]
        fourth_nucleotide=bulge[3]
        fifth_nucleotide=bulge[4]
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_potential_nt_third_nt=return_list_non_identical_nucleotides(third_nucleotide)
        list_potential_nt_fourth_nt=return_list_non_identical_nucleotides(fourth_nucleotide)
        list_potential_nt_fifth_nt=return_list_non_identical_nucleotides(fifth_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                for nucleotide_3 in list_potential_nt_third_nt:
                    for nucleotide_4 in list_potential_nt_fourth_nt:
                        for nucleotide_5 in list_potential_nt_fifth_nt:
                          new_bulge=nucleotide_1+nucleotide_2+nucleotide_3+nucleotide_4+nucleotide_5
                          list_all_potential_combinations.append(new_bulge)     
    #bulge=6nt           
    if len(bulge)==6:
        first_nucleotide=bulge[0]
        second_nucleotide=bulge[1]
        third_nucleotide=bulge[2]
        fourth_nucleotide=bulge[3]
        fifth_nucleotide=bulge[4]
        sixth_nucleotide=bulge[5]        
        list_potential_nt_first_nt=return_list_non_identical_nucleotides(first_nucleotide)
        list_potential_nt_second_nt=return_list_non_identical_nucleotides(second_nucleotide)
        list_potential_nt_third_nt=return_list_non_identical_nucleotides(third_nucleotide)
        list_potential_nt_fourth_nt=return_list_non_identical_nucleotides(fourth_nucleotide)
        list_potential_nt_fifth_nt=return_list_non_identical_nucleotides(fifth_nucleotide)
        list_potential_nt_sixth_nt=return_list_non_identical_nucleotides(sixth_nucleotide)
        list_all_potential_combinations=[]
        for nucleotide_1 in list_potential_nt_first_nt:
            for nucleotide_2 in list_potential_nt_second_nt:
                for nucleotide_3 in list_potential_nt_third_nt:
                    for nucleotide_4 in list_potential_nt_fourth_nt:
                        for nucleotide_5 in list_potential_nt_fifth_nt:
                            for nucleotide_6 in list_potential_nt_sixth_nt:
                              new_bulge=nucleotide_1+nucleotide_2+nucleotide_3+nucleotide_4+nucleotide_5+nucleotide_6
                              list_all_potential_combinations.append(new_bulge)                           
    return(list_all_potential_combinations)



        

def asemble_final_iSBHs_from_iSBHs_without_bulges_and_variable_bulge_sizes(list_potential_iSBHs_no_bulges,desired_iSBH_fold):
  if len(list_potential_iSBHs_no_bulges) >0:
    list_final_iSBH_candidates_with_bulges=[]
    list_final_folding_probabilities=[]
    
    #identifying desired bulge positions
    mis_match_intervals=identifying_corresponding_MM_positions(desired_iSBH_fold)[0]
    
    #intermediate_fold_design
    desired_intermediate_fold=determine_mfe_structure(list_potential_iSBHs_no_bulges[0], 'rna')
    
    for mismatch in mis_match_intervals:
        MM_giraffe_star=mismatch[0]
        start_interval_MM_giraffe_star=int(MM_giraffe_star[0])
        end_interval_MM_giraffe_star=int(MM_giraffe_star[-1])+1
        MM_giraffe=mismatch[1]
        start_interval_MM_giraffe=int(MM_giraffe[0])
        end_interval_MM_giraffe=int(MM_giraffe[-1])+1
              
         
        
        #re-updating iSBH structure with every MM iteration
        
        desired_intermediate_fold=desired_intermediate_fold[:start_interval_MM_giraffe_star]+('.')*len(MM_giraffe_star)+desired_intermediate_fold[end_interval_MM_giraffe_star:start_interval_MM_giraffe]+('.')*len(MM_giraffe)+desired_intermediate_fold[end_interval_MM_giraffe:]

        #print len(desired_intermediate_fold), desired_intermediate_fold

        
        list_potential_iSBHs=[] 
        for preliminary_iSBH in list_potential_iSBHs_no_bulges:

            mis_match=preliminary_iSBH[start_interval_MM_giraffe:end_interval_MM_giraffe]

            list_all_MM_combinations=return_all_bulge_combinations_except_input_for_up_to_6nt_seq(mis_match)

            list_all_bulge_combinations_for_an_intermediary_iSBH=[]
            list_all_folding_probabilities_bulge_combinations_for_an_intermediary_iSBH=[]
            for MM_element in list_all_MM_combinations:
                expanded_iSBH=preliminary_iSBH[:start_interval_MM_giraffe]+MM_element+preliminary_iSBH[end_interval_MM_giraffe:]
            
                if determine_mfe_structure(expanded_iSBH, 'rna')== desired_intermediate_fold:
                   list_all_bulge_combinations_for_an_intermediary_iSBH.append(expanded_iSBH)
                   list_all_folding_probabilities_bulge_combinations_for_an_intermediary_iSBH.append(determine_probability_secondary_structure(expanded_iSBH,desired_intermediate_fold,'rna'))
                   index_best_fold=list_all_folding_probabilities_bulge_combinations_for_an_intermediary_iSBH.index(max(list_all_folding_probabilities_bulge_combinations_for_an_intermediary_iSBH))
                   iSBH_to_keep_evolving=list_all_bulge_combinations_for_an_intermediary_iSBH[index_best_fold]
                   #print iSBH_to_keep_evolving
                   list_potential_iSBHs.append(iSBH_to_keep_evolving)
        
        list_potential_iSBHs_no_bulges=list_potential_iSBHs
       
    for final_iSBH in list_potential_iSBHs_no_bulges:
        if "TTTT" not in final_iSBH:
          list_final_iSBH_candidates_with_bulges.append(final_iSBH)
          folding_probability_final_iSBH=determine_probability_secondary_structure(final_iSBH,desired_iSBH_fold, 'rna')
          list_final_folding_probabilities.append(folding_probability_final_iSBH)
    
    return(list_final_iSBH_candidates_with_bulges, list_final_folding_probabilities)
  else:
    return('NA')

def determine_consecutive_base_pairing_interactions_between_iSBHs_no_toehold_and_trigger(iSBH_sequence, trigger_sequence, min_trigger_length):
    
    reverse_complementary_full_trigger=get_DNA_reverse_complement(trigger_sequence)
    
    seed=''
    index=1
    while index < (len(iSBH_sequence)-20):
        if (iSBH_sequence[(index-1):(index+20)] not in reverse_complementary_full_trigger) and (iSBH_sequence[(index):(index+20)] in reverse_complementary_full_trigger):
            start_nucleotides=iSBH_sequence[index:(index+20)]
            seed=seed+start_nucleotides
        if (len(seed)>=19) :
           if ((seed+ iSBH_sequence[index]) in reverse_complementary_full_trigger) and ((seed+ iSBH_sequence[index]) in iSBH_sequence):
               seed=seed+ iSBH_sequence[index]    
        index=index+1

    full_iSBH_sensing_sequence=seed
    trigger_sub_sequence=get_DNA_reverse_complement(full_iSBH_sensing_sequence)
    return(full_iSBH_sensing_sequence, len(full_iSBH_sensing_sequence), trigger_sub_sequence)
 
    
def append_toeholds_to_iSBHs_with_bulges(list_final_iSBH_candidates_with_bulges, desired_iSBH_fold_without_toehold, desired_toehold_length, desired_length_bp_tracr, trigger, folding_probability_cut_off):
   
    tracr_sequence='GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC'
    
    tracr_nt_to_bp=tracr_sequence[0:desired_length_bp_tracr]
    reverse_complement_tracr_nt_to_bp=get_DNA_reverse_complement(tracr_nt_to_bp)
    
    partial_list_iSBHs_with_bulges_and_toeholds=[]
    
    final_list_iSBHs_with_bulges_and_toeholds=[]
    list_probabilities_of_folding=[]
    list_total_ranges_of_interaction_with_trigger=[]
    
    desired_fold_iSBH_with_toehold='.'+(desired_toehold_length*'.')+(desired_length_bp_tracr*'(')+desired_iSBH_fold_without_toehold+(desired_length_bp_tracr*')')
    
    for iSBH_seq in list_final_iSBH_candidates_with_bulges:
      trigger_sub_sequence=determine_consecutive_base_pairing_interactions_between_iSBHs_no_toehold_and_trigger(iSBH_seq, trigger, 24)[2]
      ranges_of_interaction_with_trigger=len(trigger_sub_sequence)+ desired_toehold_length
      list_total_ranges_of_interaction_with_trigger.append(ranges_of_interaction_with_trigger)
      
      index=0
      list_toehold_donor_sequences=[]
      list_length_toehold_donor_sequences=[]
      while index <=len(trigger):
        potential_sub_sequence=trigger[index:(index+len(trigger_sub_sequence))]   
        if potential_sub_sequence == trigger_sub_sequence:
            toehold_design_space=trigger[(index+len(trigger_sub_sequence)):]
            list_toehold_donor_sequences.append(toehold_design_space)
            list_length_toehold_donor_sequences.append(len(toehold_design_space))
        index=index+1
    
      index_toehold_donor_sequence=list_length_toehold_donor_sequences.index(max(list_length_toehold_donor_sequences))
      final_toehold_donnor_sequence=get_DNA_reverse_complement(list_toehold_donor_sequences[index_toehold_donor_sequence])
    
      index_2=0
      while index_2 <= (len(final_toehold_donnor_sequence)-desired_toehold_length):
        potential_toehold=final_toehold_donnor_sequence[index_2:(index_2+desired_toehold_length)]
        
        potential_iSBH_with_toehold='G'+potential_toehold+reverse_complement_tracr_nt_to_bp+str(iSBH_seq)+tracr_nt_to_bp
        
        if determine_mfe_structure(potential_iSBH_with_toehold, 'rna')==desired_fold_iSBH_with_toehold:
            partial_list_iSBHs_with_bulges_and_toeholds.append(potential_iSBH_with_toehold)
        
        index_2=index_2+1
    
      for toeholded_iSBH in partial_list_iSBHs_with_bulges_and_toeholds:
          folding_probability=determine_probability_secondary_structure(toeholded_iSBH, desired_fold_iSBH_with_toehold,'rna')
          if folding_probability > folding_probability_cut_off:
              final_list_iSBHs_with_bulges_and_toeholds.append(toeholded_iSBH)
              list_probabilities_of_folding.append(folding_probability)
     
    return(final_list_iSBHs_with_bulges_and_toeholds, list_probabilities_of_folding, list_total_ranges_of_interaction_with_trigger, desired_fold_iSBH_with_toehold)

def selecting_no_toehold_iSBHs_that_still_maintain_final_folding_when_adding_tracr(list_iSBHs_no_toehold_no_GC, desired_iSBH_fold_no_GC):
    list_validated_iSBHs=[]
    tracr_extra_fold=".......((((....)))).((((.......))))............((((....)))).((((((...))))))"
    tracr_sequence="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
    desired_final_fold=".("+desired_iSBH_fold_no_GC+')'+tracr_extra_fold
    for iSBH_sequence in list_iSBHs_no_toehold_no_GC:
      full_iSBH_sequence='GC'+iSBH_sequence+tracr_sequence
      if determine_mfe_structure(full_iSBH_sequence, 'rna')==desired_final_fold:
          list_validated_iSBHs.append(iSBH_sequence)
    
    return (list_validated_iSBHs)

def selecting_iSBHs_with_toehold_that_still_maintain_final_folding_when_adding_tracr(list_iSBHs_with_toehold_and_GC, desired_iSBH_fold_no_tracr):
    list_validated_iSBHs=[]
    tracr_extra_fold=".......((((....)))).((((.......))))............((((....)))).((((((...))))))"
    tracr_sequence="GTTTTAGAGCTAGAAATAGCAAGTTAAAATAAGGCTAGTCCGTTATCAACTTGAAAAAGTGGCACCGAGTCGGTGC"
    
    desired_final_fold=desired_iSBH_fold_no_tracr+tracr_extra_fold
    for iSBH_sequence in list_iSBHs_with_toehold_and_GC:
      full_iSBH_sequence=iSBH_sequence+tracr_sequence[1:]
      full_iSBH_measured_fold= determine_mfe_structure(full_iSBH_sequence, 'rna')
      if full_iSBH_measured_fold ==desired_final_fold:
          list_validated_iSBHs.append(iSBH_sequence)
    
    return (list_validated_iSBHs)


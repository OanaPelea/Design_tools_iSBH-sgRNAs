# -*- coding: utf-8 -*-
"""
Created on Tue May 01 14:50:08 2018

@author: Pelea Oana
"""

from iSBH_functions import ISBH_split
from iSBH_functions import reverse_complement
from iSBH_functions import generate_plasmid_map
from iSBH_functions import standard_sequence_format

# Input the name of the file where you provided full-length iSBH sequences (FW strands)
input_file_name='SL14_giraffe_designs.txt'

output_1_name="to_order_"+ input_file_name.replace('.txt','.csv')
output_2_name="Records_oligoes_"+input_file_name
output_3_name='Seq_records_'+input_file_name.replace('.txt','.csv')

# Provide the name of a desired output file
input_file_folder='inputs/input_iSBH_sequences/'
output_files_CSV_folder='outputs/oligos_to_order/'
output_files_records_folder='outputs/oligos_record_files/'
output_folder_plasmid_maps='outputs/plasmid_maps/'

input_file= input_file_folder+ input_file_name
to_order_output_file=open(output_files_CSV_folder + output_1_name, "w")
record_output_file=open(output_files_records_folder + output_2_name, "w")
sequencing_summary=open(('outputs/sequencing_tables/'+output_3_name), 'w')

#Input all design characteristics

length_spacer=int(raw_input("Length spacer: "))
length_Giraffe_extension=int(raw_input("Length giraffe extension: "))
length_iSBH_loop=int(raw_input("Length sensing loop: "))
length_Giraffe_extension_star=int(raw_input("Length Giraffe extension * : "))
length_spacer_star=int(raw_input("Length spacer * :"))
total_iSBH_length=length_spacer+length_Giraffe_extension+length_iSBH_loop+length_Giraffe_extension_star+length_spacer_star
additional_C_beginning_spacer_star="YES" #YES/NO
trigger_length=length_spacer+length_Giraffe_extension+length_iSBH_loop
CTS_PAM_site="AGG"
my_name_initials='OP_'
experiment_name=''
output_split_iSBH_sequences=str(raw_input("Output split iSBH sequences? #YES/NO : ")) #"YES"/"NO"
output_spacer_sequences=str(raw_input("Output spacer sequences? #YES/NO : ")) #"YES"/"NO"
output_triger_sequences=str(raw_input("Output trigger sequences? #YES/NO : ")) #"YES"/"NO"
output_CTS_sequences=str(raw_input("Output split 1xCTS sequences? #YES/NO : ")) #"YES"/"NO"

#Additional inputs that could change manually                                  
                                   
extra_G_U6_expression="YES" #"YES"/"NO"
output_plasmid_maps="YES" #"YES"/"NO"

#
if additional_C_beginning_spacer_star == 'YES':
  my_extra_C="C"
  my_extra_G="G"
else:
  my_extra_C=""
  my_extra_G=""  

if extra_G_U6_expression== "YES":
  U6_G="G"
  U6_C="C"
else:
  U6_G=""
  U6_C=""

index_2=int(raw_input("Index of the iSBH sequence column in the input file: "))
  
#Checking if all sequences provided respect design specifications
index=int(raw_input("Index first plasmid: "))
for line in open(input_file).readlines():    
    f=line.split()
    iSBH_name=f[0]
    iSBH_label=my_name_initials+experiment_name+iSBH_name    
    iSBH_sequence=f[index_2]
    print iSBH_sequence
    
    if len(iSBH_sequence) != total_iSBH_length:
        print " "
        print "For", iSBH_name,"double-check design specifications"
        
    else:
        print " "
        print "For", iSBH_name, "provided design specifications are right"      
        complementary_iSBH_right_orientation=reverse_complement(iSBH_sequence)
        
        
        #defining important sequences 

        my_spacer_star=iSBH_sequence[0:length_spacer_star]
        my_Giraffe_extension_star=iSBH_sequence[length_spacer_star:(length_spacer_star+length_Giraffe_extension_star)]
        my_loop=iSBH_sequence[(length_spacer_star+length_Giraffe_extension_star):(length_spacer_star+length_Giraffe_extension_star+length_iSBH_loop)]
        my_Giraffe_extension=iSBH_sequence[(length_spacer_star+length_Giraffe_extension_star+length_iSBH_loop):(length_spacer_star+length_Giraffe_extension_star+length_iSBH_loop+length_Giraffe_extension)]
        my_spacer=iSBH_sequence[(length_spacer_star+length_Giraffe_extension_star+length_iSBH_loop+length_Giraffe_extension):(length_spacer_star+length_Giraffe_extension_star+length_iSBH_loop+length_Giraffe_extension+length_spacer)]
        first_nucleotide_spacer_star=my_spacer_star[0]
        first_nucleotide_spacer=my_spacer[0]
        
   
        #designing iSBH FW oligos       
        if first_nucleotide_spacer_star == "G":
          if additional_C_beginning_spacer_star == 'YES':
            Bbs_5_prime_FW_strand = "CACC"+U6_G+my_extra_C
          else:
            Bbs_5_prime_FW_strand = "CACC"+my_extra_C
        else:
          Bbs_5_prime_FW_strand = "CACC"+U6_G+my_extra_C  
        full_iSBH_seq_FW=Bbs_5_prime_FW_strand+iSBH_sequence
        

        #designing iSBH RW oligos               
        if first_nucleotide_spacer_star == "G":
            if additional_C_beginning_spacer_star == 'YES':
              term_3_prime_RW_strand=my_extra_G+U6_C
            else:
              term_3_prime_RW_strand=my_extra_G
        else:
            term_3_prime_RW_strand=my_extra_G+U6_C
        BBs_5_prime_RW_strand="AAAC"
        full_iSBH_seq_RW=BBs_5_prime_RW_strand+complementary_iSBH_right_orientation+term_3_prime_RW_strand

        first_fw_oligo_1=ISBH_split(full_iSBH_seq_FW, full_iSBH_seq_RW)[0]
        first_rw_oligo_2=ISBH_split(full_iSBH_seq_FW, full_iSBH_seq_RW)[1]
        second_fw_oligo_3=ISBH_split(full_iSBH_seq_FW, full_iSBH_seq_RW)[2]      
        second_rw_oligo_4=ISBH_split(full_iSBH_seq_FW, full_iSBH_seq_RW)[3]
        
        #generating plasmid map iSBHs
        if (output_plasmid_maps=="YES") and (output_split_iSBH_sequences=="YES"):
            plasmid_map_name=output_folder_plasmid_maps+str(index)+'_'+iSBH_label+'.dna'
            plasmid_map_file=open(plasmid_map_name, 'w')
            
            top_oligo=standard_sequence_format(first_fw_oligo_1)+standard_sequence_format(second_fw_oligo_3)
            bottom_oligo=standard_sequence_format(second_rw_oligo_4)+standard_sequence_format(first_rw_oligo_2)
            plasmid_map_file.write(generate_plasmid_map('U6_sgRNA', 'BbsI', 'BbsI', top_oligo, bottom_oligo))
            plasmid_map_file.close()
            sequencing_summary.write(plasmid_map_name)
            sequencing_summary.write(',')
            sequencing_summary.write(top_oligo)
            sequencing_summary.write('\n')
            index=index+1    
        #designing spacer FW
        if first_nucleotide_spacer == "G": 
          spacer_fw="CACC"+my_spacer
        else:
          spacer_fw="CACC"+U6_G+my_spacer     

        
        #designing spacer RW
        reverse_complement_spacer=reverse_complement(my_spacer)
        
        if first_nucleotide_spacer == "G":
          spacer_rw="AAAC"+reverse_complement_spacer
        else:
          spacer_rw="AAAC"+reverse_complement_spacer+U6_C
        
        #generating plasmid map spacer
        if (output_plasmid_maps=="YES") and (output_spacer_sequences=="YES"):
            plasmid_map_name_spa=output_folder_plasmid_maps+str(index)+'_'+iSBH_label+'_SPA.dna'
            plasmid_map_spa_file=open(plasmid_map_name_spa, 'w')
            
            plasmid_map_spa_file.write(generate_plasmid_map('U6_sgRNA', 'BbsI', 'BbsI', spacer_fw, spacer_rw))
            plasmid_map_spa_file.close()
            
            sequencing_summary.write(plasmid_map_name_spa)
            sequencing_summary.write(',')
            sequencing_summary.write(spacer_fw)
            sequencing_summary.write('\n')
            
            index=index+1  
        #designing trigger RW
        trig_rw=my_spacer_star+my_Giraffe_extension_star+my_loop
        trig_fw=reverse_complement(trig_rw)  
        full_trig_rw="GCAT"+trig_rw
        full_trig_fw="CCAC"+trig_fw
        
        #generating plasmid map trigger
        if (output_plasmid_maps=="YES") and (output_triger_sequences=="YES"):
            plasmid_map_name_trig=output_folder_plasmid_maps+str(index)+'_'+iSBH_label+'_TRIG.dna'
            plasmid_map_trig_file=open(plasmid_map_name_trig, 'w')
            plasmid_map_trig_file.write(generate_plasmid_map('U6_TM2', 'BbsI', 'BbsI', full_trig_fw, full_trig_rw))
            plasmid_map_trig_file.close()
            
            sequencing_summary.write(plasmid_map_name_trig)
            sequencing_summary.write(',')
            sequencing_summary.write(full_trig_fw)
            sequencing_summary.write('\n')

            index=index+1

        #designing 1xCTS FW
        CTS=my_spacer+CTS_PAM_site        
        CTS_fw="CTAGA"+my_spacer+CTS_PAM_site+"GG"
       
        #designing 1XCTS RW     
        CTS_rw="CGCGCC"+reverse_complement(CTS)+"T"

        #generating plasmid map 1xCTS
        if (output_plasmid_maps=="YES") and (output_CTS_sequences=="YES"):
            plasmid_map_name_CTS=output_folder_plasmid_maps+str(index)+'_'+iSBH_label+'_1xCTS.dna'
            plasmid_map_CTS_file=open(plasmid_map_name_CTS, 'w')
            plasmid_map_CTS_file.write(generate_plasmid_map('1XCTS_CFP', 'XbaI', 'AscI', CTS_fw, CTS_rw))
            plasmid_map_CTS_file.close()
            
            sequencing_summary.write(plasmid_map_name_CTS)
            sequencing_summary.write(',')
            sequencing_summary.write(CTS_fw)
            sequencing_summary.write('\n')
            index=index+1
        
        #printing in the to_order file
        if output_split_iSBH_sequences == "YES":

          to_order_output_file.write(str(iSBH_label+"_o1 "))
          to_order_output_file.write(',')
          to_order_output_file.write(first_fw_oligo_1)
          to_order_output_file.write('\n')
          to_order_output_file.write(str(iSBH_label+"_o2 "))
          to_order_output_file.write(',')
          to_order_output_file.write(first_rw_oligo_2)
          to_order_output_file.write('\n')
          to_order_output_file.write(str(iSBH_label+"_o3 "))
          to_order_output_file.write(',')
          to_order_output_file.write(second_fw_oligo_3)
          to_order_output_file.write('\n')       
          to_order_output_file.write(str(iSBH_label+"_o4 "))
          to_order_output_file.write(',')
          to_order_output_file.write(second_rw_oligo_4)
          to_order_output_file.write('\n')
          
        if output_spacer_sequences == "YES":
          to_order_output_file.write(str(iSBH_label+"_spac_fw "))
          to_order_output_file.write(',')
          to_order_output_file.write(spacer_fw)
          to_order_output_file.write('\n')
          to_order_output_file.write(str(iSBH_label+"_spac_rw "))
          to_order_output_file.write(',')
          to_order_output_file.write(spacer_rw)
          to_order_output_file.write('\n') 
        
        if output_triger_sequences == "YES":
          to_order_output_file.write(str(iSBH_label+"_trig_fw"))
          to_order_output_file.write(',')
          to_order_output_file.write(full_trig_fw)
          to_order_output_file.write('\n')       
          to_order_output_file.write(str(iSBH_label+"_trig_rw"))
          to_order_output_file.write(',')
          to_order_output_file.write(full_trig_rw)
          to_order_output_file.write('\n')
          
        if output_CTS_sequences == "YES":
          to_order_output_file.write(str(iSBH_label+"_CTS_fw"))
          
          to_order_output_file.write(',')
          to_order_output_file.write(CTS_fw)
          to_order_output_file.write('\n')
          to_order_output_file.write(str(iSBH_label+"_CTS_rw"))
          to_order_output_file.write(',')
          to_order_output_file.write(CTS_rw)
          to_order_output_file.write('\n')
          to_order_output_file.write('\n')
          print str(iSBH_label+"_CTS_fw"), CTS_fw
          print str(iSBH_label+"_CTS_rw"), CTS_rw
          
        
        #printing in the output_record        
        record_output_file.write(iSBH_name)
        record_output_file.write('_full_length_iSBH_FW (no restriction sites): ')
        record_output_file.write(my_extra_C+iSBH_sequence)
        record_output_file.write( "\n")
        
        record_output_file.write(iSBH_name)
        record_output_file.write('_spacer_FW (no restriction sites): ')
        record_output_file.write(my_spacer)
        record_output_file.write( "\n")
        
        record_output_file.write(iSBH_name)
        record_output_file.write('_trigger_fw(no restriction sites):')
        record_output_file.write(trig_fw)
        record_output_file.write( "\n")
        
        record_output_file.write(iSBH_name)
        record_output_file.write('_CTS_fw(no restriction sites):')
        record_output_file.write(CTS)
        record_output_file.write( "\n")
        record_output_file.write( "\n")
        record_output_file.write( "\n")
        
        
record_output_file.write( "\n")
record_output_file.write( "\n")
record_output_file.write("##################################################################################################################")
record_output_file.write( "\n")
record_output_file.write( "\n")        
record_output_file.write( "The following sequences have been designed according to your input specifications:")
record_output_file.write( "\n")
record_output_file.write( "     Total iSBH length: ")
record_output_file.write(str(total_iSBH_length))
record_output_file.write( "\n")
record_output_file.write( "     Spacer length: ")
record_output_file.write( str(length_spacer))
record_output_file.write( "\n")
record_output_file.write("     Giraffe extension length: ")
record_output_file.write(str(length_Giraffe_extension))
record_output_file.write( "\n")
record_output_file.write( "     Loop length:")
record_output_file.write(str(length_iSBH_loop))
record_output_file.write( "\n")
record_output_file.write( "     Giraffe extension* length: ")
record_output_file.write(str(length_Giraffe_extension_star))
record_output_file.write( "\n")
record_output_file.write( "     Spacer* length: ")
record_output_file.write( str(length_spacer_star))
record_output_file.write( "\n")
if additional_C_beginning_spacer_star == 'YES':
  record_output_file.write( "     You have an additional C at the beginning of your spacer* sequence")
  record_output_file.write( "\n")
if extra_G_U6_expression== "YES":
  record_output_file.write( "     For iSBHs, spacers (BUT NO triggers!!!) that do not start with G, an extra G has been added at the beginning of the sequence to facilitate transcription from U6 promoter")
  record_output_file.write( "\n")
else:
  record_output_file.write( "    An extra G for improving U6 transcription has NOT been added. \n")
record_output_file.write(  "     Trigger length:")
record_output_file.write( str(trigger_length))
record_output_file.write( "\n")
record_output_file.write("     Selected PAM site for the 1X CTS sequence: ")
record_output_file.write(CTS_PAM_site)
record_output_file.write( "\n")
record_output_file.write( "     iSBH oligos are ready to be cloned in pcDNA3.1-U6_sgRNA_6xT-SV40_iBue_PA (BbsI restriction sites)")
record_output_file.write( "\n")
record_output_file.write( "     spacer oligos are ready to be cloned in pcDNA3.1-U6_sgRNA_6xT-SV40_iBue_PA (BbsI restriction sites)")
record_output_file.write( "\n")
record_output_file.write( "     trigger oligos are ready to be cloned in U6-TM2Emp-6T_iBlue (BbsI restriction sites)")
record_output_file.write( "\n")
record_output_file.write( "     CTS oligos are ready to be cloned in p035_pause-HBG-CFP-pA (XbaI and AscI restriction sites)")
record_output_file.write( "\n")
record_output_file.write("##################################################################################################################")
record_output_file.write( "\n")
record_output_file.write( "\n")
        
record_output_file.close()
to_order_output_file.close()
sequencing_summary.close()
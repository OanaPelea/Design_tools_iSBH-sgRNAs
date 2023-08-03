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

input_filename='inputs/input_spacers.csv'
output_filename='outputs/optimal_spacer_folds.csv'

output_file=open(output_filename,'w')

output_file.write('sgRNA_name , spacer_sequence , spacer_star_sequence , minimal_fold_sequence,,,,\n')
list_spacers=[] 
for line in open (input_filename).readlines():
    f=line.split(',')
    sgRNA_name=f[0]
    sgRNA_sequence=standard_sequence_format(f[1])
    list_spacers.append(sgRNA_sequence)

print list_spacers

#determining most suitable spacer * sequences
list_spacer_stars_and_spacers=[]

index=0
for spacer in list_spacers:
    index=index+1
    spacer_star=determine_most_stable_backfold_for_given_spacer(spacer)[0]
    SPA=determine_most_stable_backfold_for_given_spacer(spacer)[1]
    small_iSBH=determine_most_stable_backfold_for_given_spacer(spacer)[2]
    sgRNA_name='sgRNA_'+str(index)
    print sgRNA_name, SPA, spacer_star, small_iSBH
    output_line=str(sgRNA_name+','+ SPA+','+ spacer_star+','+small_iSBH +','+ ',,,\n')
    output_file.write(output_line)
output_file.close()
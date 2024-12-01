import numpy as np
import pandas as pd
import os
import json
import time

from utilities import ncbi_data
from utilities_polycistron import (mirna_polycistron_data, 
                                   mirna_members_data,
                                   check_motif,
                                   mirna_loading)


# Load config

with open(os.path.join('config.json'), 'r') as f:
    config = json.load(f)

with open(os.path.join('files', 'input_data.json'), 'r') as f:
    config_input_data = json.load(f)

compl_dict = config['compl_dict']
file_name = config_input_data['general_information']['file_name']
folder_name = config_input_data['general_information']['output_folder_name']
output_folder =  os.path.join('files', 'outputs', folder_name)

# Load information from user config
ncbi_name = config_input_data['gene_information']['ncbi_name']
cluster_name = config_input_data['polycistron_data']['cluster_name']
mirna_members = config_input_data['polycistron_data']['mirna_members']
pair_mode = config_input_data['polycistron_data']['pair_mode']

# Load information from settings config
mirna_polycistrons_file = config['polycistron_data']['mirna_polycistrons_file']
mirna_in_polycistron_file = config['polycistron_data']['mirna_in_polycistron_file']
rna_fold_file = config['polycistron_data']['rna_fold_file']
rna_fold_out_file = config['polycistron_data']['rna_fold_out_file']
rna_cofold_file = config['polycistron_data']['rna_cofold_file']
rna_cofold_out_file = config['polycistron_data']['rna_cofold_out_file']

# File with selected sequences
mirna_file = os.path.join(output_folder, f'{file_name}_mirna_with_features.csv')


# Files for control results
if not os.path.exists(os.path.join(output_folder, 'sirna_with_structure_error.txt')):
    with open(os.path.join(output_folder, 'sirna_with_structure_error.txt'), 'w') as f:
        f.write('')
if not os.path.exists(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt')):
    with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt'), 'w') as f:
        f.write('')

# Get transcript_sequence
refseq_sequence, _, _, _, _ = ncbi_data(ncbi_name)

# Get natural polycistron data
(members, 
 scaffold_mirna, 
 left_flank, 
 right_flank, 
 sequence) = mirna_polycistron_data(mirna_polycistrons_file, 
                                    cluster_name, 
                                    flank_size = 200)

# Get matural mirna
scaffold, all_members, mirna_in_polycistrons = mirna_members_data(mirna_in_polycistron_file, 
                                                                  mirna_members, 
                                                                  scaffold_mirna, 
                                                                  members)

# Check neccesary parameters and save
loops_motifs_check, spacer_motifs_check = check_motif(scaffold_mirna, sequence, all_members)

motif_results = pd.concat([loops_motifs_check, spacer_motifs_check])
motif_results.to_csv(os.path.join(output_folder, f'{cluster_name}_motif_check_results.csv'), index=None)

# Load information for mirna from mirbase db
mirna_with_info = mirna_loading(mirna_in_polycistrons)


time.sleep(100000)


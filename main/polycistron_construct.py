"""Code for polycistron processing"""

import pandas as pd
import os
import json
import time
from datetime import date

from utilities import ncbi_data
from utilities_polycistron import (mirna_polycistron_data, 
                                   mirna_members_data,
                                   check_motif,
                                   mirna_loading,
                                   sirna_data,
                                   mirna_sirna_pairs,
                                   get_new_scaffold,
                                   full_scaffold_structures,
                                   hairpin_plots,
                                   save_results)
from snapgene_utilities import gene_bank_file, find_elements_polycistron

DATE_TODAY = str(date.today())

"""Load config"""

with open(os.path.join('config.json'), 'r') as f:
    config = json.load(f)

with open(os.path.join('files', 'input_data.json'), 'r') as f:
    config_input_data = json.load(f)

compl_dict = config['compl_dict']
file_name = config_input_data['general_information']['file_name']
folder_name = config_input_data['general_information']['output_folder_name']
output_folder =  os.path.join('files', 'outputs', folder_name)
vienna_output_directory = os.path.join('files', 'outputs', 'vienna')

"""Load information from user config"""
ncbi_name = config_input_data['gene_information']['ncbi_name']
cluster_name = config_input_data['polycistron_data']['cluster_name']
mirna_members = config_input_data['polycistron_data']['mirna_members']
pair_mode = config_input_data['polycistron_data']['pair_mode']

"""Load information from settings config"""
mirna_polycistrons_file = config['polycistron_data']['mirna_polycistrons_file']
mirna_in_polycistron_file = config['polycistron_data']['mirna_in_polycistron_file']

"""Data for Vienna RNA tool"""
rna_fold_file = os.path.join(vienna_output_directory, 'RNAfold_data.fa')
rna_fold_out_file = os.path.join(vienna_output_directory, 'RNAfold_output.out')
rna_cofold_file = os.path.join(vienna_output_directory, 'data_for_RNAcofold.fasta')
rna_cofold_out_file = os.path.join(vienna_output_directory, 'RNAcofoldresults.txt')

"""File with selected sequences"""
mirna_file = os.path.join(output_folder, f'{file_name}_mirna_with_features.csv')

"""Files for control results"""
if not os.path.exists(os.path.join(output_folder, 'sirna_with_structure_error.txt')):
    with open(os.path.join(output_folder, 'sirna_with_structure_error.txt'), 'w') as f:
        f.write('')

with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt'), 'w') as f:
    f.write('')

"""Get transcript_sequence"""
refseq_sequence, _, _, _, _ = ncbi_data(ncbi_name)

"""Get natural polycistron data"""
(members, 
 scaffold_mirna, 
 left_flank, 
 right_flank, 
 sequence) = mirna_polycistron_data(mirna_polycistrons_file, 
                                    cluster_name, 
                                    flank_size = 200)

"""Get matural mirna"""
scaffold, all_members, mirna_in_polycistrons = mirna_members_data(mirna_in_polycistron_file, 
                                                                  mirna_members, 
                                                                  scaffold_mirna, 
                                                                  members)

"""Check neccesary parameters and save"""
loops_motifs_check, spacer_motifs_check = check_motif(scaffold_mirna, sequence, all_members)

motif_results = pd.concat([loops_motifs_check, spacer_motifs_check])
motif_results.to_csv(os.path.join(output_folder, f'{cluster_name}_motif_check_results.csv'), index=None)

"""Load information for mirna from mirbase db"""
mirna_with_info = mirna_loading(mirna_in_polycistrons)

"""Load new mirna/sirna data"""
all_sequences = sirna_data(mirna_file, refseq_sequence, compl_dict, output_folder)

"""Compare sequences of natural and new mirna/sirna"""
pairs, all_alignments = mirna_sirna_pairs(mirna_in_polycistrons, scaffold, all_sequences, type=pair_mode, k=0.4)

"""Compute folding with new scaffold"""
(scaffold_clean, scaffold_new, 
 all_old_sequences, all_new_sequences, 
 mirna_names, sirna_names, 
 all_old_structure, all_new_structure) = get_new_scaffold(scaffold, mirna_in_polycistrons, 
                                                        all_alignments, pairs, 
                                                        rna_fold_file, rna_fold_out_file, 
                                                        rna_cofold_file, rna_cofold_out_file,
                                                        output_folder,
                                                        vienna_output_directory,
                                                        compl_dict)

new_scaffold_seq = left_flank + scaffold_new + right_flank
with open(os.path.join(output_folder, f'{file_name}_polycistron_sequence.fa'), 'w') as f:
    f.write(f'> {file_name}\n')
    f.write(new_scaffold_seq + '\n')

"""Compute all new sequences"""
full_scaffold_structures(scaffold_clean, scaffold_new, left_flank, right_flank, 
                         all_old_sequences, all_new_sequences, 
                         mirna_names, sirna_names, sequence,
                         rna_fold_file, rna_fold_out_file,
                        vienna_output_directory, output_folder)

"""Make hairpin pictures"""
hairpin_plots(all_old_sequences, 
              all_new_sequences, 
              mirna_names, 
              sirna_names, 
              rna_fold_file, 
              rna_fold_out_file,
              vienna_output_directory,
              output_folder)

"""Save all sequences"""
save_results(output_folder, 
            all_old_structure, 
            all_new_structure, 
            mirna_names, 
            sirna_names,
            all_sequences,
            all_new_sequences,
            ncbi_name,
            refseq_sequence)

"""Make Gene Bank file for SnapGene tool"""
elements_list, oligos = find_elements_polycistron(new_scaffold_seq, 
                                                    all_sequences, 
                                                    sirna_names, 
                                                    all_new_sequences)

gene_bank_file(ncbi_name, 
               new_scaffold_seq, 
               DATE_TODAY, 
               elements_list, 
               f'{file_name}_polycistron', 
               output_folder,
               oligos=oligos)
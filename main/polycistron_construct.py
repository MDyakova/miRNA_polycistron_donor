import numpy as np
import pandas as pd
import os
import json,
import time

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






time.sleep(100000)


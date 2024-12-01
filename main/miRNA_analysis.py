import numpy as np
import pandas as pd
import os
import json 
import subprocess
import time
from utilities import (ensembl_data, 
                       ncbi_data, 
                       data_for_vienna, 
                       rnacofold_results, 
                       rnafold_results, 
                       rnafold_mirna_results, 
                       correct_sequences,
                       transcript_features,
                       compute_GC_context)

# Load config

with open(os.path.join('config.json'), 'r') as f:
    config = json.load(f)

with open(os.path.join('files', 'input_data.json'), 'r') as f:
    config_input_data = json.load(f)

compl_dict = config['compl_dict']
file_name = config_input_data['general_information']['file_name']
folder_name = config_input_data['general_information']['output_folder_name']

output_folder =  os.path.join('files', 'outputs', folder_name)
os.makedirs(output_folder, exist_ok=True)

# Load Gene ID, transcript ID and cromosome name from config

gene_information = config_input_data['gene_information']

gene_name = gene_information['gene_name']
ncbi_name = gene_information['ncbi_name']
chromosome_name = gene_information['chromosome_name']

# Get gene coordinates 
start_gene, end_gene = ensembl_data(gene_name)

# Get transript features
refseq_sequence, features_list, cds_start, cds_end, exons = ncbi_data(ncbi_name)

# Load possible miRNA (excel table)
mirna_sequences_path = config_input_data['data']['mirna_table']
mirna_data = pd.read_excel(mirna_sequences_path)
mirna_data = mirna_data[pd.notna(mirna_data['sequence'])]

# miRNA sequences preparation (compute number of relevant sources for each miRNA)
mirna_data['sequence'] = mirna_data['sequence'].apply(lambda p: p.upper().replace('U', 'T'))
mirna_data = mirna_data.fillna('')
mirna_data['count_sourses'] = mirna_data['sourse']

def sourse_agg(sourse_data):
    return ' '.join(sourse_data)

mirna_data = mirna_data.groupby(by=['sequence'], as_index=False).agg({'Name/gene name':max, 'sourse': sourse_agg, 
                                                        'count_sourses':len, 'type':max, 'Note':max})
mirna_data = mirna_data[['Name/gene name', 'sequence', 'sourse', 'count_sourses', 'type', 'Note']]

# Check that in table correct strand and fix if not (miRNA should be complimentary to transcript)
mirna_correct = correct_sequences(mirna_data, refseq_sequence, compl_dict)

# Save correct sequences to file
mirna_correct_sequences_path = os.path.join(output_folder, 
                                            f'{file_name}_correct_sequences.csv')
mirna_correct.to_csv(mirna_correct_sequences_path, index=None)

# Directory for all vienna results
vienna_output_directory = os.path.join('files', 'outputs', 'vienna')
os.makedirs(vienna_output_directory, exist_ok=True)

# Make files for Vienna RNA cofold
(rnacofold_input, 
 alignments_path,
 transcript_path,
 mirna_seq_path) = data_for_vienna(mirna_correct, 
                                    refseq_sequence, 
                                    file_name, 
                                    ncbi_name)

# Launch RNAcofold tool
rnacofold_data = rnacofold_results(rnacofold_input, vienna_output_directory)

# Launch RNAfold
all_alignments_results = pd.read_csv(alignments_path)
rnafold_data = rnafold_results(mirna_correct, 
                               all_alignments_results, 
                               transcript_path, 
                               vienna_output_directory, 
                               ncbi_name)

# Launch RNAfold for mirna
mirnafold_data = rnafold_mirna_results(mirna_seq_path, vienna_output_directory)

# Add transcript features
mirna_with_features = transcript_features(mirna_correct, cds_start, cds_end, exons)

# Join all results together
mirna_with_features['GC_content'] = mirna_with_features['sequence'].apply(lambda p:compute_GC_context(p))
mirna_with_features = pd.merge(mirna_with_features, rnacofold_data, on=['sequence'], how='left')
mirna_with_features = pd.merge(mirna_with_features, rnafold_data, on=['sequence'], how='left')
mirna_with_features = pd.merge(mirna_with_features, mirnafold_data, on=['sequence'], how='left')
mirna_with_features['choice'] = 0

# Save results
mirna_with_features.to_csv(os.path.join(output_folder,
                                        f'{file_name}_mirna_with_features.csv'), 
                                        index=None)

time.sleep(100000)






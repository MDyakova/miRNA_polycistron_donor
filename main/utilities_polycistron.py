import numpy as np
import pandas as pd
import os
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from Bio import Entrez
from Bio import SeqIO
import matplotlib.pyplot as plt
import subprocess

def mirna_polycistron_data(file_name, cluster_name, flank_size = 200):
    with open(file_name) as f:
        mirna_polycistrons = f.read()
    
    mirna_polycistrons = mirna_polycistrons.split('>')
    
    sequence = list(filter(lambda p: cluster_name in p, mirna_polycistrons))[0]
    members = sequence.split('\n')[0].split(' ')[2:]
    sequence = sequence.split('\n')[1]
    left_flank = sequence.split('[')[0][-flank_size:]
    right_flank = sequence.split(']')[-1][:flank_size]
    scaffold =  sequence.split(left_flank)[1].split(right_flank)[0]
    return members, scaffold, left_flank, right_flank, sequence

def mirna_members_data(mirna_in_polycistron_file, mirna_members, scaffold, members):
    mirna_in_polycistrons = pd.read_csv(mirna_in_polycistron_file)
    
    mirna_in_polycistrons = mirna_in_polycistrons[mirna_in_polycistrons['mir_id'].apply(lambda p: p in members)]
    mirna_in_polycistrons = mirna_in_polycistrons.sort_values(by=['start'])
    mirna_in_polycistrons.reset_index(inplace=True)
    
    all_members = list(mirna_in_polycistrons['mirbase_name'])
    
    mirna_in_polycistrons['in_list'] = mirna_in_polycistrons['mirbase_name'].apply(lambda p: 1 if p in mirna_members else 0)
    mirna_in_polycistrons['index'] = mirna_in_polycistrons.index
    scaffold_copy = scaffold
    
    mirna_not_in_polycistron = mirna_in_polycistrons[mirna_in_polycistrons['in_list']==0]
    for ind in mirna_not_in_polycistron['index']:
        hairpin_i = scaffold_copy.split('[')[ind+1].split(']')[0]
        hairpin_i = scaffold_copy.split(hairpin_i)[0][-1:] + hairpin_i + scaffold.split(hairpin_i)[1][:1]
        scaffold = scaffold.replace(hairpin_i, '')
    
    mirna_in_polycistrons = mirna_in_polycistrons[mirna_in_polycistrons['in_list']==1]

    return scaffold, all_members, mirna_in_polycistrons

def check_motif(scaffold, sequence, all_members):

    loop_motifs = ['TGT', 'GTG']
    
    loops = [seq.split(']')[0].split('}')[1].split('{')[0] for seq in scaffold.split('[')[1:]]
    
    loops_motifs_check = []
    for loop_seq in loops:
        is_correct = False
        for motif in loop_motifs:
            if motif in loop_seq:
               is_correct = True 
        loops_motifs_check.append([loop_seq, is_correct])
    
    loops_motifs_check = pd.DataFrame(loops_motifs_check, columns=('sequence', 'with_motif'))
    loops_motifs_check['mirna'] = all_members
    loops_motifs_check['type'] = 'loop motif'
    # mirna without loop motifs
    loops_motifs_check = loops_motifs_check[loops_motifs_check['with_motif']==False]
    
    spacers = [seq.split('[')[0][:40] for seq in sequence.split(']')[1:]]
    
    spacer_motifs_check = []
    
    for spacer in spacers:
        is_correct = np.max([((spacer[i]=='C') & (spacer[i+3]== 'C'))for i in range(len(spacer)-3)])
        spacer_motifs_check.append([spacer, is_correct])
        # break
    
    spacer_motifs_check = pd.DataFrame(spacer_motifs_check, columns=('sequence', 'with_motif'))
    spacer_motifs_check['mirna'] = all_members
    spacer_motifs_check['type'] = 'spacer motif'
    # mirna without spacer motifs
    spacer_motifs_check = spacer_motifs_check[spacer_motifs_check['with_motif']==False]

    # pd.concat([loops_motifs_check, spacer_motifs_check]).to_csv(output_folder + cluster_name + '_motif_check_results.csv', index=None)
    return loops_motifs_check, spacer_motifs_check

def source_number(evidence):
    sources = 0
    for string in evidence.split('[')[1:]:
        for string_i in string.split(']')[0].split(','):
            if string_i.isdigit():
               sources += 1
            else:
                sources += int(string_i.split('-')[1]) - int(string_i.split('-')[0]) + 1
    return sources

def mirna_loading(mirna_in_polycistrons):
    mirna_in_polycistrons['evidence1_n'] = mirna_in_polycistrons['evidence1'].apply(lambda p: source_number(p))
    mirna_in_polycistrons['evidence2_n'] = mirna_in_polycistrons['evidence2'].apply(lambda p: source_number(p))
    
    mirna_in_polycistrons['mirna_strand'] = mirna_in_polycistrons['evidence1_n'] - mirna_in_polycistrons['evidence2_n']
    mirna_in_polycistrons['mirna_strand'] = mirna_in_polycistrons['mirna_strand'].apply(lambda p: 0 if p>=0 else 1)

    return mirna_in_polycistrons

def sirna_data(sirna_file, refseq_sequence, comp_dict, output_folder):

    with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt')) as f:
        sirna_with_incorrect_energy = f.read()
    sirna_with_incorrect_energy = sirna_with_incorrect_energy.split('\n')
    
    if '.csv' in sirna_file:
        sirnas = pd.read_csv(sirna_file)
    else:
        sirnas = pd.read_excel(sirna_file)
    sirnas = sirnas[sirnas['choice']==1]
    
    if 'name' not in sirnas.columns:
        sirnas['name'] = 'si' + sirnas['Name/gene name'] + '_' + sirnas['start_mirna'].astype('str')
        
    sirnas = sirnas[sirnas['name'].apply(lambda p: p not in sirna_with_incorrect_energy)]

    # sirnas = sirnas.head(5)
    
    all_sequences = []
    for i in range(len(sirnas)):
        seq = sirnas.iloc[i]['sequence']
        rev_seq = ''.join([comp_dict[s] for s in seq][::-1])
    
        all_aligh1 = pairwise2.align.localms(refseq_sequence, seq, 2, -3, -5, -2)
        score1 = all_aligh1[0].score
        
        all_aligh2 = pairwise2.align.localms(refseq_sequence, rev_seq, 2, -3, -5, -2)
        score2 = all_aligh2[0].score
    
        if score1>=score2:
            sirna_seq = rev_seq
            all_aligh = all_aligh1[0]
        else:
            sirna_seq = seq
            all_aligh = all_aligh2[0]
    
        
        # start = all_aligh.start
        # end = all_aligh.end
        delta = 23 - (all_aligh.end - all_aligh.start)
        gene_seq = all_aligh.seqA[all_aligh.start - delta:all_aligh.end]
        sirna_seq = ''.join([comp_dict[s] for s in gene_seq][::-1])
        all_sequences.append([sirnas['name'].iloc[i], sirna_seq])
    
    all_sequences = []
    for i in range(len(sirnas)):
        seq = sirnas.iloc[i]['sequence']
        rev_seq = ''.join([comp_dict[s] for s in seq][::-1])
    
        all_aligh1 = pairwise2.align.localms(refseq_sequence, seq, 2, -3, -5, -2)
        score1 = all_aligh1[0].score
        
        all_aligh2 = pairwise2.align.localms(refseq_sequence, rev_seq, 2, -3, -5, -2)
        score2 = all_aligh2[0].score
    
        if score1>=score2:
            sirna_seq = rev_seq
            all_aligh = all_aligh1[0]
        else:
            sirna_seq = seq
            all_aligh = all_aligh2[0]
        
        # start = all_aligh.start
        # end = all_aligh.end
        delta = 23 - (all_aligh.end - all_aligh.start)
        gene_seq = all_aligh.seqA[all_aligh.start - delta:all_aligh.end]
        sirna_seq = ''.join([comp_dict[s] for s in gene_seq][::-1])
        all_sequences.append([sirnas['name'].iloc[i], sirna_seq])

    return all_sequences

def mirna_sirna_pairs(mirna_in_polycistrons, scaffold, all_sequences, type = 'similar', k=0.4):
    all_alignments = []
    for i in range(np.minimum(len(mirna_in_polycistrons), len(all_sequences))):
        mirna_name = mirna_in_polycistrons.iloc[i]['mirbase_name']
        strand = mirna_in_polycistrons.iloc[i]['mirna_strand']
        mirna_seq = scaffold.split('[')[1:][i].split('{')[1:][strand].split('}')[0]
        hairpin_seq = scaffold.split('[')[1:][i].split(']')[0].replace('{', '').replace('}', '')
        for sirna_name, sirna_seq in all_sequences:
            delta = len(mirna_seq) - len(sirna_seq)
            all_alignments.append([mirna_name, sirna_name, mirna_seq, sirna_seq, hairpin_seq,
                                   np.sum([(s1 == s2) for s1, s2 in zip(mirna_seq[:len(mirna_seq)-delta], sirna_seq)])])
    all_alignments = pd.DataFrame(all_alignments, columns=('mirna_name', 'sirna_name', 'mirna_seq', 'sirna_seq', 'hairpin_seq', 'score'))
    all_alignments_g = all_alignments.groupby(by=['sirna_name'], as_index=False).max()[['sirna_name', 'score']]
    all_alignments = pd.merge(all_alignments, all_alignments_g, on=['sirna_name'], suffixes=('', '_max'))
    all_alignments['k'] = all_alignments['score']/all_alignments['score_max']
    all_alignments_f = all_alignments[all_alignments['k']>=k]
    all_alignments_c = all_alignments_f.groupby(by=['mirna_name'], as_index=False).count()[['sirna_name', 'mirna_name']]
    all_alignments_c = all_alignments_c.sort_values(by=['sirna_name'])
        
    if type == 'similar':
        sirna_pairs = []
        mirna_pairs = []
        for mirna_name in all_alignments_c['mirna_name']:
            all_alignments_i = all_alignments[(all_alignments['mirna_name']==mirna_name) 
                                                & (all_alignments['sirna_name'].apply(lambda p: p not in sirna_pairs))]
            all_alignments_i = all_alignments_i.sort_values(by=['score'], ascending=False)
            sirna_name = all_alignments_i.iloc[0]['sirna_name']
            mirna_pairs.append(mirna_name)
            sirna_pairs.append(sirna_name)
        pairs = pd.DataFrame(data={'sirna_names':sirna_pairs, 'mirna_names':mirna_pairs})
    else:
        sirna_pairs = []
        mirna_pairs = []
        for mirna_name in all_alignments_c['mirna_name']:
            all_alignments_i = all_alignments_f[(all_alignments_f['mirna_name']==mirna_name) 
                                                & (all_alignments_f['sirna_name'].apply(lambda p: p not in sirna_pairs))]
            all_alignments_i = all_alignments_i.sample(frac=1)
            sirna_name = all_alignments_i.iloc[0]['sirna_name']
            mirna_pairs.append(mirna_name)
            sirna_pairs.append(sirna_name)
        pairs = pd.DataFrame(data={'sirna_names':sirna_pairs, 'mirna_names':mirna_pairs})
    return pairs, all_alignments

def new_sequence_function(mirna_structure_list, sirna_seq, compl_dict, delta = 1):
    new_sequence = []
    mirna_checked = 0
    sirna_pos = 0
    for step in range(len(mirna_structure_list[1])):
        if (mirna_structure_list[1][step].islower()) | (mirna_structure_list[0][step].islower()):
            new_sequence.append([mirna_structure_list[0][step], 
                         mirna_structure_list[1][step], 
                         mirna_structure_list[2][step], 
                         mirna_structure_list[3][step],
                         mirna_structure_list[4][step]])
        elif (mirna_structure_list[1][step].isupper()) | (mirna_structure_list[0][step].isupper()):
            if mirna_checked < delta:
                new_sequence.append([mirna_structure_list[0][step], 
                             mirna_structure_list[1][step], 
                             mirna_structure_list[2][step], 
                             mirna_structure_list[3][step],
                             mirna_structure_list[4][step]])
            elif  mirna_checked >= (len(sirna_seq) + delta):
                new_sequence.append([mirna_structure_list[0][step], 
                             mirna_structure_list[1][step], 
                             mirna_structure_list[2][step], 
                             mirna_structure_list[3][step],
                             mirna_structure_list[4][step]])
            else:
                if mirna_structure_list[2][step]=='|':
                    new_sequence.append([mirna_structure_list[0][step], 
                                         sirna_seq[sirna_pos], 
                                         mirna_structure_list[2][step], 
                                         compl_dict[sirna_seq[sirna_pos]],
                                         mirna_structure_list[4][step]])
                    sirna_pos += 1
                else:
                    if mirna_structure_list[4][step]=='-':
                        new_sequence.append([sirna_seq[sirna_pos], 
                                             mirna_structure_list[1][step], 
                                             mirna_structure_list[2][step], 
                                             mirna_structure_list[3][step],
                                             '-'])
                        sirna_pos += 1
                    else:
                        new_sequence.append([sirna_seq[sirna_pos], 
                                             mirna_structure_list[1][step], 
                                             mirna_structure_list[2][step], 
                                             mirna_structure_list[3][step],
                                             sirna_seq[sirna_pos]])
                        sirna_pos += 1
                
            mirna_checked += 1
                
        else:
            new_sequence.append([mirna_structure_list[0][step], 
                         mirna_structure_list[1][step], 
                         mirna_structure_list[2][step], 
                         mirna_structure_list[3][step],
                         mirna_structure_list[4][step]])
    new_sequence_list = []
    for seq in np.array(new_sequence).T:
        new_sequence_list.append(''.join(seq))
    return new_sequence, new_sequence_list

def list_to_seq(mirna_structure_list):
    new_sequence = ''
    for step in range(len(mirna_structure_list[1])):
        if mirna_structure_list[0][step] < mirna_structure_list[1][step]:
            new_sequence += mirna_structure_list[1][step]
        else:
            new_sequence += mirna_structure_list[0][step]
    
    new_sequence += mirna_structure_list[2][step]

    for step in range(len(mirna_structure_list[1])-1, -1, -1):
        if mirna_structure_list[3][step] < mirna_structure_list[4][step]:
            new_sequence += mirna_structure_list[4][step]
        else:
            new_sequence += mirna_structure_list[3][step]
    
    new_sequence_r = new_sequence.replace('-', '').upper().replace('U', 'T')
    new_sequence = new_sequence.upper().replace('U', 'T')
    return new_sequence, new_sequence_r

def letter_size(mirna_structure, new_sequence_list):
    mirna_structure_split = mirna_structure.split('\n')
    new_sequence_list_split = [list(i) for i in new_sequence_list]
    if len(mirna_structure_split[-1]) < len(mirna_structure_split[1]):
        mirna_structure_split[-1] = mirna_structure_split[-1] + ' '
    for step in range(len(new_sequence_list[1])):
        if mirna_structure_split[0][step].islower():
            new_sequence_list_split[0][step] = new_sequence_list_split[0][step].lower()
        if mirna_structure_split[1][step].islower():
            new_sequence_list_split[1][step] = new_sequence_list_split[1][step].lower()
        if mirna_structure_split[2][step].islower():
            new_sequence_list_split[2][step] = new_sequence_list_split[2][step].lower()
        if mirna_structure_split[3][step].islower():
            new_sequence_list_split[3][step] = new_sequence_list_split[3][step].lower()
        if mirna_structure_split[4][step].islower():
            new_sequence_list_split[4][step] = new_sequence_list_split[4][step].lower()
    return [''.join(i) for i in new_sequence_list_split]

def get_new_sequence(mirna_structure_list, 
                     mirna_structure, 
                     sirna_seq, 
                     hairpin_seq, 
                     strand, 
                     rna_fold_file, 
                     rna_fold_out_file,
                     vienna_output_directory, 
                     compl_dict,
                     delta = 1):
    if strand == 0:
            new_sequence, new_sequence_list = new_sequence_function(mirna_structure_list, sirna_seq, compl_dict, delta = delta)
    if strand == 1:
        mirna_structure_list = [i[::-1] for i in mirna_structure_list][::-1]
        new_sequence, new_sequence_list = new_sequence_function(mirna_structure_list, sirna_seq, compl_dict, delta = delta)
        new_sequence_list = [i[::-1] for i in new_sequence_list][::-1]
        # break

   
    new_sequence_list = letter_size(mirna_structure, new_sequence_list)
    new_sequence, new_sequence_r = list_to_seq(new_sequence_list)

    with open(rna_fold_file, 'w') as f:
        f.write('> old seq \n')
        f.write(hairpin_seq.replace('T', 'U') + '\n')
        f.write('> new seq \n')
        f.write(new_sequence_r.replace('T', 'U')  + '\n')

    # !RNAfold -p -d2 < data/RNAfold_data.fa > data/RNAfold_output.out
    # Launch RNAfold
    with open(rna_fold_file, "r") as infile, open(rna_fold_out_file, "w") as outfile:
        subprocess.run(["RNAfold", "-p", "-d2"], stdin=infile, stdout=outfile, check=True, cwd=vienna_output_directory)

    with open(rna_fold_out_file) as f:
        RNAfold_output = f.readlines()
    old_structure = RNAfold_output[2].split(' ')[0]
    new_structure = RNAfold_output[8].split(' ')[0]
    
    old_energy = float(RNAfold_output[2].split(' ')[1].split('(')[1].split(')')[0])
    new_energy = float(RNAfold_output[8].split(' ')[1].split('(')[1].split(')')[0])
    
    delta_energy = (old_energy - new_energy)
    
    delta_structures = old_structure == new_structure

    structure_compare = [int(old_structure[i] == new_structure[i]) for i in range(len(old_structure))]
    structure_compare = np.diff(np.where(np.array(structure_compare)==1)[0]) - 1
    
    return new_sequence, new_sequence_list, structure_compare, delta_energy

def get_mature_sequence_coords(sequence):
    start_nt0_list = ([i for i in range(len(list(sequence[0]))) if list(sequence[0])[i].isupper()])
    start_nt1_list = ([i for i in range(len(list(sequence[1]))) if list(sequence[1])[i].isupper()])
    if len(start_nt0_list)>0:
        start_nt0 = np.min(start_nt0_list)
    else:
        # value more that length of mirna
        start_nt0 = 100
    if len(start_nt1_list)>0:
        start_nt1 = np.min(start_nt1_list)
    else:
        # value more that length of mirna
        start_nt1 = 100

    end_nt0_list = ([i for i in range(len(list(sequence[3]))) if list(sequence[3])[i].isupper()])
    end_nt1_list = ([i for i in range(len(list(sequence[4]))) if list(sequence[4])[i].isupper()])
    if len(end_nt0_list)>0:
        end_nt0 = np.max(end_nt0_list)
    else:
        # value more that length of mirna
        end_nt0 = -1
        
    if len(end_nt1_list)>0:
        end_nt1 = np.max(end_nt1_list)
    else:
        # value less that 0
        end_nt1 = -1
          
    start_nt = np.minimum(start_nt0, start_nt1)
    end_nt = np.maximum(end_nt0, end_nt1) + 1

    return start_nt, end_nt

def energy_sites(mature_seq_list, strand, rna_cofold_file, rna_cofold_out_file, vienna_output_directory, size = 7):
    if strand == 0:
        left_side = [i[:size] for i in mature_seq_list]
        right_side = [i[-size:] for i in mature_seq_list]
    else:
        right_side = [i[:size] for i in mature_seq_list]
        left_side = [i[-size:] for i in mature_seq_list]
    
    left_seq1 = ''.join([list(filter(lambda p: p!=' ', [a, b]))[0] for a,b in zip(left_side[0], left_side[1])])
    left_seq2 = ''.join([list(filter(lambda p: p!=' ', [a, b]))[0] for a,b in zip(left_side[3], left_side[4])])[::-1]
    
    right_seq1 = ''.join([list(filter(lambda p: p!=' ', [a, b]))[0] for a,b in zip(right_side[0], right_side[1])])
    right_seq2 = ''.join([list(filter(lambda p: p!=' ', [a, b]))[0] for a,b in zip(right_side[3], right_side[4])])[::-1]
    
    with open(rna_cofold_file, 'w') as f:
        f.write('> sample left \n')
        f.write(left_seq1 + '&' + left_seq2 + '\n')
        f.write('> sample right \n')
        f.write(right_seq1 + '&' + right_seq2 + '\n')
    
    # !RNAcofold -p < $rna_cofold_file > $rna_cofold_out_file
    with open(rna_cofold_file, "r") as infile, open(rna_cofold_out_file, "w") as outfile:
        subprocess.run(["RNAcofold", "-p"], stdin=infile, stdout=outfile, check=True, cwd=vienna_output_directory)
    
    with open(rna_cofold_out_file) as f:
        RNAcofoldresults = f.readlines()
    
    left_energy = float(RNAcofoldresults[4].split('delta G binding=')[-1].split('\n')[0])
    right_energy = float(RNAcofoldresults[9].split('delta G binding=')[-1].split('\n')[0])

    return left_energy, right_energy

def get_new_scaffold(scaffold, 
                     mirna_in_polycistrons, 
                     all_alignments, pairs, 
                     rna_fold_file, 
                     rna_fold_out_file, 
                     rna_cofold_file, 
                     rna_cofold_out_file,
                     output_folder,
                     vienna_output_directory,
                     compl_dict):
    scaffold_clean = scaffold.replace('{', '').replace('}', '')
    scaffold_new = scaffold.replace('{', '').replace('}', '')
    all_old_sequences = []
    all_new_sequences = []
    mirna_names = []
    sirna_names = []
    all_old_structure = []
    all_new_structure = []

    with open(os.path.join(output_folder, 'sirna_with_structure_error.txt')) as f:
        sirna_structure_e = f.read()
    with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt')) as f:
        sirna_incorrect_e = f.readlines()
    for mirna_name in pairs['mirna_names']:
        mirna_data = mirna_in_polycistrons[mirna_in_polycistrons['mirbase_name']==mirna_name]
        mirna_structure = mirna_data['structure'].max()
        mirna_seq = all_alignments[all_alignments['mirna_name']==mirna_name]['mirna_seq'].max()
        hairpin_seq = all_alignments[all_alignments['mirna_name']==mirna_name]['hairpin_seq'].max()
        strand = mirna_in_polycistrons[mirna_in_polycistrons['mirbase_name']==mirna_name]['mirna_strand'].max()
    
        sirna_name = pairs[pairs['mirna_names']==mirna_name]['sirna_names'].max()
        sirna_seq = all_alignments[all_alignments['sirna_name']==sirna_name]['sirna_seq'].max()
    
        delta = len(mirna_seq) - len(sirna_seq)
        mirna_structure_list = mirna_structure.split('\n')
        if len(mirna_structure_list[-1]) < len(mirna_structure_list[1]):
            mirna_structure_list[-1] = mirna_structure_list[-1] + ' '
    
        (new_sequence, 
         new_sequence_list, 
         structure_compare, 
         delta_energy) = get_new_sequence(mirna_structure_list, mirna_structure, 
                                        sirna_seq, hairpin_seq, strand, 
                                        rna_fold_file,  
                                        rna_fold_out_file,
                                        vienna_output_directory, compl_dict, delta = 1)
        
        if (np.max(structure_compare) >= 5) | (len(np.where(structure_compare>0)[0])>3) | (delta_energy<-5):
            print('check another delta '  +  mirna_name)
    
            (new_sequence, 
             new_sequence_list, 
             structure_compare, 
             delta_energy) = get_new_sequence(mirna_structure_list, mirna_structure, 
                                                sirna_seq, hairpin_seq, strand, 
                                                rna_fold_file, 
                                                rna_fold_out_file, 
                                                vienna_output_directory, compl_dict, delta = 0)
            
            if (np.max(structure_compare) >= 5) | (len(np.where(structure_compare>2)[0])>3) | (delta_energy<-5):
                with open(os.path.join(output_folder, 'sirna_with_structure_error.txt')) as f:
                    sirna_structure_e = f.read()
                if (sirna_name + ' ' + mirna_name) not in sirna_incorrect_e:
                    with open(os.path.join(output_folder, 'sirna_with_structure_error.txt'), 'a') as f:
                        f.write(sirna_name + ' ' + mirna_name + '\n')
                # print('Structure error: ' + sirna_name + ' ' + mirna_name)

        sequence = new_sequence_list.copy()
        start_nt, end_nt = get_mature_sequence_coords(sequence)
        
        mature_seq_list = [i[start_nt:end_nt] for i in sequence]
        
        left_energy, right_energy = energy_sites(mature_seq_list, 
                                                 strand, 
                                                 rna_cofold_file, 
                                                 rna_cofold_out_file, 
                                                 vienna_output_directory)
        
        
        if right_energy - left_energy > 0:
            with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt')) as f:
                sirna_incorrect_e = f.readlines()
            if sirna_name not in sirna_incorrect_e:
                with open(os.path.join(output_folder, 'sirna_with_incorrect_energy.txt'), 'a') as f:
                    f.write(sirna_name + '\n')
            # print(sirna_name + ' sites have incorrect energies')
    
        if len(scaffold_new.split(hairpin_seq))==2:
            scaffold_new = scaffold_new.replace(hairpin_seq, new_sequence.replace('-', ''))
        # else:
        #     'Scaffold sequence error ' + mirna_name
        # break
    
        all_old_sequences.append(hairpin_seq)
        all_new_sequences.append(new_sequence.replace('-', ''))
        mirna_names.append(mirna_name)
        sirna_names.append(sirna_name)
        all_old_structure.append(mirna_structure.split('\n'))
        all_new_structure.append(new_sequence_list)
    
    scaffold_clean = scaffold_clean.replace('[', '').replace(']', '')
    scaffold_new = scaffold_new.replace('[', '').replace(']', '')
    return scaffold_clean, scaffold_new, all_old_sequences, all_new_sequences, mirna_names, sirna_names, all_old_structure, all_new_structure

def structure_prob(file_name):

    with open(file_name) as f:
        dp_data = f.read()
    
    probabilities = pd.DataFrame([i.split(' ') for i in dp_data.split('%start of base pair probability data\n')[1].split('\nshowpage')[0].split('\n')], 
                                                columns=('pos_1', 'pos_2', 'prob', 'type'))
    
    probabilities['prob'] = probabilities['prob'].astype('float')
    probabilities['pos_1'] = probabilities['pos_1'].astype('int')
    probabilities['pos_2'] = probabilities['pos_2'].astype('int')
    
    probabilities_1 = probabilities[probabilities['type']=='ubox'].groupby(by=['pos_1'], as_index=False).max()
    probabilities_1['pos'] = probabilities_1['pos_1']
    
    probabilities_2 = probabilities[probabilities['type']=='ubox'].groupby(by=['pos_2'], as_index=False).max()
    probabilities_2['pos'] = probabilities_2['pos_2']
    
    probabilities = pd.concat([probabilities_1, probabilities_2])
    probabilities = probabilities.groupby(by=['pos'], as_index=False).max().sort_values('pos')
    
    sequence = ''.join(dp_data.split('/sequence')[1].split(') }')[0].split('\\\n')[1:-1])
    
    probabilities_all = pd.DataFrame(data={'pos':np.arange(1, len(sequence)+1)})
    probabilities_all = pd.merge(probabilities_all, probabilities, on=['pos'], how='left').fillna(0)
    
    structure = ''.join(probabilities_all['prob'].apply(lambda p: '|' if p>0.85 else '.'))

    return structure, probabilities_all, sequence

def structure_coords(file_name, probabilities_all):
    
    with open(file_name) as f:
        ss_data = f.read()

    coords = ss_data.split('/coor [\n')[1].split('\n] def')[0].split('\n')
    coords = [[float(c[1:-1].split(' ')[0]), float(c[1:-1].split(' ')[1])] for c in coords]

    probabilities_all['X'] = np.array(coords).T[0]
    probabilities_all['Y'] = np.array(coords).T[1]

    return probabilities_all

def structure_picture(probabilities_all, all_sequences, all_names, sequence, picture_name):
    colors = plt.cm.jet(np.linspace(0,1,len(all_sequences)))
    plt.figure()
    plt.plot(probabilities_all['X'], probabilities_all['Y'], 'k', linewidth=0.5)
    for step, (seq, name) in enumerate(zip(all_sequences, all_names)):
        seq_start = len(sequence.split(seq.replace('T', 'U'))[0])
        seq_end = seq_start + len(seq)
    
        prob_seq = probabilities_all[(probabilities_all['pos']>=seq_start) & (probabilities_all['pos']<seq_end)]
    
        plt.plot(prob_seq['X'], prob_seq['Y'], c = colors[step], linewidth=1, label=name)
        plt.title(picture_name)
        plt.legend()
        plt.savefig(output_folder + output_name + ' ' + picture_name + '.jpg')

def full_scaffold_structures(scaffold_clean, scaffold_new, left_flank, right_flank, 
                             all_old_sequences, all_new_sequences,  
                             mirna_names, sirna_names, sequence):
    scaffold_clean = left_flank + scaffold_clean + right_flank
    scaffold_new = left_flank + scaffold_new + right_flank
    
    with open(rna_fold_file, 'w') as f:
        f.write('> old seq \n')
        f.write(scaffold_clean.replace('T', 'U') + '\n')
        f.write('> new seq \n')
        f.write(scaffold_new.replace('T', 'U')  + '\n')
    
    # !RNAfold -p -d2 < $rna_fold_file > $rna_fold_out_file
    
    with open(rna_fold_out_file) as f:
        RNAfold_output = f.readlines()
    old_structure = RNAfold_output[2].split(' ')[0]
    new_structure = RNAfold_output[8].split(' ')[0]
    
    old_energy = float(RNAfold_output[2].split(' ')[1].split('(')[1].split(')')[0])
    new_energy = float(RNAfold_output[8].split(' ')[1].split('(')[1].split(')')[0])
    
    file_name = 'old_dp.ps'
    structure, probabilities_all, sequence = structure_prob(file_name)
    
    file_name = 'old_ss.ps'
    probabilities_all = structure_coords(file_name, probabilities_all)

    structure_picture(probabilities_all, all_old_sequences, mirna_names, sequence, 'mirna polycistron folding')

    file_name = 'new_dp.ps'
    structure, probabilities_all, sequence = structure_prob(file_name)
    
    file_name = 'new_ss.ps'
    probabilities_all = structure_coords(file_name, probabilities_all)

    structure_picture(probabilities_all, all_new_sequences, sirna_names, sequence, 'sirna polycistron folding')


def hairpin_plots(all_old_sequences, all_new_sequences, mirna_names, sirna_names):
    colors = plt.cm.jet(np.linspace(0,1,len(mirna_names)))
    fig, axs = plt.subplots(len(mirna_names)*2, 1, figsize=(5,3 * len(mirna_names)))
    pos = 0
    for step, (seq_old, seq_new, mirna_name, sirna_name) in enumerate(zip(all_old_sequences, all_new_sequences, mirna_names, sirna_names)):  
        with open(rna_fold_file, 'w') as f:
            f.write('> old seq \n')
            f.write(seq_old.replace('T', 'U') + '\n')
            f.write('> new seq \n')
            f.write(seq_new.replace('T', 'U')  + '\n')
        
        # !RNAfold -p -d2 < data/RNAfold_data.fa > data/RNAfold_output.out
        
        with open(rna_fold_out_file) as f:
            RNAfold_output = f.readlines()
        old_structure = RNAfold_output[2].split(' ')[0]
        new_structure = RNAfold_output[8].split(' ')[0]
        
        old_energy = float(RNAfold_output[2].split(' ')[1].split('(')[1].split(')')[0])
        new_energy = float(RNAfold_output[8].split(' ')[1].split('(')[1].split(')')[0])
    
        file_name = 'old_dp.ps'
        structure, probabilities_all, sequence = structure_prob(file_name)
        
        file_name = 'old_ss.ps'
        probabilities_all = structure_coords(file_name, probabilities_all)

        axs[pos].plot(probabilities_all['Y'], probabilities_all['X'], 'k', linewidth=0.5)
        axs[pos].title.set_text(mirna_name + ' energy = ' + str(old_energy))
        pos += 1
        
        file_name = 'new_dp.ps'
        structure, probabilities_all, sequence = structure_prob(file_name)
        
        file_name = 'new_ss.ps'
        probabilities_all = structure_coords(file_name, probabilities_all)
    
        axs[pos].plot(probabilities_all['Y'], probabilities_all['X'], c = colors[step], linewidth=1)
        axs[pos].title.set_text(sirna_name + ' energy = ' + str(new_energy))
        pos += 1
    
    fig.tight_layout(pad=1.0)
    plt.savefig(output_folder + output_name + ' ' + 'pri-mirna folding' + '.jpg')


"""Functions for main scripts"""

import numpy as np
import pandas as pd
import os
import json
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from pyensembl import EnsemblRelease
import seaborn as sns
from Bio import Entrez
from Bio import SeqIO
import subprocess

# Load config

with open(os.path.join('config.json'), 'r') as f:
    config = json.load(f)

compl_dict = config['compl_dict']

def ensembl_data(gene_name):
    data = EnsemblRelease(109)
    gene = data.genes_by_name(gene_name)[0]
    start_gene = gene.start - 200
    end_gene = gene.end + 200
    return start_gene, end_gene

def ncbi_data(ncbi_name):
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="nucleotide", id=ncbi_name, rettype="gb", retmode="text")
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)

    features_list = []
    exons = []
    exon_name = 1
    for feature in seq_record.features:
        features_list.append([feature.type, int(feature.location.start+1), int(feature.location.end), 
                        feature.qualifiers])
        if feature.type == 'CDS':
            cds_start = int(feature.location.start+1)
            cds_end = int(feature.location.end)
        if feature.type == 'exon':
            exons.append(['exon' + str(exon_name), int(feature.location.start+1), int(feature.location.end)])
            exon_name += 1

    return refseq_sequence, features_list, cds_start, cds_end, exons

def data_for_vienna(mirna_data, refseq_sequence, file_name, ncbi_name):
    all_alignments = []
    for step, sirna_i in enumerate(mirna_data['sequence']):
        sirna_i = sirna_i.upper().replace('U', 'T')
        sirna_r_i = ''.join([compl_dict[s] for s in sirna_i][::-1])
        if (sirna_i[:-2] in refseq_sequence):
            all_aligh = pairwise2.align.localms(refseq_sequence, sirna_i, 2, -3, -5, -2)

            sirna_start = np.argmax([l!='-' for l in all_aligh[0].seqB])
            all_alignments.append([sirna_i, sirna_r_i, refseq_sequence[sirna_start:sirna_start + len(sirna_i)],
                                ])

        elif (sirna_r_i[2:] in refseq_sequence):
            all_aligh_r = pairwise2.align.localms(refseq_sequence, sirna_r_i, 2, -3, -5, -2)

            sirna_start = np.argmax([l!='-' for l in all_aligh_r[0].seqB])
            all_alignments.append([sirna_i, sirna_i, refseq_sequence[sirna_start:sirna_start + len(sirna_i)],
                                ])

        else:
            all_aligh = pairwise2.align.localms(refseq_sequence, sirna_i, 2, -3, -5, -2)
            all_aligh_r = pairwise2.align.localms(refseq_sequence, sirna_r_i, 2, -3, -5, -2)
            if all_aligh[0].score >= all_aligh_r[0].score:
                sirna_start = np.argmax([l!='-' for l in all_aligh[0].seqB])
                all_alignments.append([sirna_i, sirna_r_i, refseq_sequence[sirna_start:sirna_start + len(sirna_i)],
                                    ])
            else:
                sirna_start = np.argmax([l!='-' for l in all_aligh_r
                                        [0].seqB])
                all_alignments.append([sirna_i, sirna_i, refseq_sequence[sirna_start:sirna_start + len(sirna_i)],
                                    ])

    all_alignments_results = pd.DataFrame(all_alignments, columns=('sirna', 'sirna_s', 'gene_seq'))
    alignments_path = os.path.join('files', 'outputs', 'vienna', f'{file_name}_for_rnacofold.csv')
    all_alignments_results.to_csv(alignments_path, index=False)

    sequences_path = os.path.join('files', 'outputs', 'vienna', f'{file_name}_for_rnacofold.fasta')
    with open(sequences_path, 'w') as f:
        for step, pair in enumerate(all_alignments):
            f.write('> ' + pair[0] + '\n')
            f.write(pair[1] + '&' + pair[2] + '\n')

    transcript_path = os.path.join('files', 'outputs', 'vienna', f'{file_name}_transcript.fasta')
    with open(transcript_path, 'w') as f:
        f.write('> ' + ncbi_name + '\n')
        f.write(refseq_sequence + '\n')

    mirna_seq_path = os.path.join('files', 'outputs', 'vienna', f'{file_name}_for_rnafold_mirna.fasta')
    with open(mirna_seq_path, 'w') as f:
        for step, pair in enumerate(all_alignments):
            f.write('> ' + pair[0] + '\n')
            f.write(pair[1] + '\n')

    return sequences_path, alignments_path, transcript_path, mirna_seq_path


def rnacofold_results(rnacofold_input, vienna_output_directory):
    # Launch RNAcofold tool
    rnacofold_output = rnacofold_input.replace('.fasta', '.txt')
    with open(rnacofold_input, "r") as infile, open(rnacofold_output, "w") as outfile:
        subprocess.run(["RNAcofold", "-p"], stdin=infile, stdout=outfile, check=True, cwd=vienna_output_directory)

    with open(rnacofold_output) as f:
        rnacofoldresults = f.readlines()

    rnacofold_data = []
    for i in range(0, len(rnacofoldresults), 5):
        rnacofold_data.append([rnacofoldresults[i].split(' ')[1].split('\n')[0], 
                            rnacofoldresults[i+4].split('delta G binding=')[1].split('\n')[0]])

    rnacofold_data = pd.DataFrame(rnacofold_data, columns=('sequence', 'delta_G_cofold'))
    rnacofold_data.drop_duplicates(subset='sequence', inplace=True)

    return rnacofold_data

def rnafold_results(mirna_data, all_alignments_results, transcript_path, vienna_output_directory, ncbi_name):
    # Launch RNAfold
    rnafold_output = transcript_path.replace('.fasta', '.out')
    with open(transcript_path, "r") as infile, open(rnafold_output, "w") as outfile:
        subprocess.run(["RNAfold", "-p", "-d2"], stdin=infile, stdout=outfile, check=True, cwd=vienna_output_directory)

    with open(os.path.join(vienna_output_directory, f'{ncbi_name}_dp.ps')) as f:
        rnafold_mrna = ''.join(f.readlines())

    mrna_seq = rnafold_mrna.split("sequence { (\\\n")[1].split('\\\n) }')[0].replace('\\\n', '').replace('U', 'T')
    bpp = rnafold_mrna.split('%start of base pair probability data\n')[1].split('\nshowpage')[0]
    bbp = pd.DataFrame([i.split(' ') for i in bpp.split('\n')], columns=('p1', 'p2', 'prob', 'type'))

    bbp = bbp[bbp['type'] == 'ubox']
    bbp['prob'] = bbp['prob'].astype('float')

    bbp['p1'] = bbp['p1'].astype('int') - 1
    bbp['p2'] = bbp['p2'].astype('int') - 1

    rnafold_data = []
    for sirna_i in mirna_data['sequence']:
        gene_s = all_alignments_results[all_alignments_results['sirna']==sirna_i]['gene_seq'].max()
        start_sirna = len(mrna_seq.split(gene_s)[0])
        end_sirna = start_sirna + len(sirna_i)
        p1 = (bbp[((bbp['p1']>=start_sirna) & (bbp['p1']<=end_sirna))])
        p2 = (bbp[((bbp['p2']>=start_sirna) & (bbp['p2']<=end_sirna))])
        pair_prob = pd.concat([p1, p2])['prob'].mean()
        rnafold_data.append([sirna_i, pair_prob])

    rnafold_data = pd.DataFrame(rnafold_data, columns=('sequence', 'pair_probability'))
    rnafold_data.drop_duplicates(subset='sequence', inplace=True)
    return rnafold_data

def rnafold_mirna_results(mirna_seq_path, vienna_output_directory):
    # Launch RNAcofold tool
    rnafold_output = mirna_seq_path.replace('.fasta', '.out')
    with open(mirna_seq_path, "r") as infile, open(rnafold_output, "w") as outfile:
        subprocess.run(["RNAfold", "-p", "-d2"], stdin=infile, stdout=outfile, check=True, cwd=vienna_output_directory)

    with open(rnafold_output) as f:
        mirnafold = f.readlines()

    mirnafold_data = []
    for i in range(0, len(mirnafold), 6):
        mirnafold_data.append([mirnafold[i].split(' ')[1].replace('U', 'T').split('\n')[0], 
                            float(mirnafold[i+2].split(' ')[-1].split(')')[0].replace('(', ''))])
        
    mirnafold_data = pd.DataFrame(mirnafold_data, columns=('sequence', 'delta_G_fold'))
    mirnafold_data.drop_duplicates(subset='sequence', inplace=True)

    return mirnafold_data

def rnahybrid_coding_results(transcripts_path, lncrna_path, mirna_seq_path, vienna_output_directory):
    # Launch RNAcofold tool
    rnahybrid_output = mirna_seq_path.replace('.fasta', '.txt')
    with open(rnahybrid_output, "w") as outfile:
        subprocess.run(
            ["RNAhybrid", "-t", transcripts_path, "-q", lncrna_path, "-s", "3utr_human"],
            stdout=outfile,
            check=True,
            cwd=vienna_output_directory
            )

    with open(rnahybrid_output) as f:
        rna_hybrid_crna = f.read()

    rna_hybrid_crna = rna_hybrid_crna.split('target:')[1:]
    rna_hybrid_crna = list(filter(lambda p: 'target too long' not in p, rna_hybrid_crna))

    rna_hybrid_crna_data = []
    for rna_i in rna_hybrid_crna:
        gene_name = rna_i.split(' ')[1].split('\n')[0]
        sirna = rna_i.split('miRNA : ')[1].split('\n')[0]
        mfe = rna_i.split('mfe: ')[1].split(' ')[0]
        p_value = rna_i.split('p-value: ')[1].split('\n')[0]
        rna_hybrid_crna_data.append([gene_name, sirna, float(mfe), float(p_value)])

    rna_hybrid_crna_data = pd.DataFrame(rna_hybrid_crna_data, 
                                        columns=('gene_name', 'sequence', 'mfe_transcript', 'p_value'))
    rna_hybrid_crna_data = rna_hybrid_crna_data[rna_hybrid_crna_data['p_value']<=0.05]
    rna_hybrid_crna_data = rna_hybrid_crna_data[rna_hybrid_crna_data['mfe_transcript']<=-25]

    rna_hybrid_crna_data = rna_hybrid_crna_data.groupby(by=['sequence'], as_index=False).agg({'mfe_transcript': 'mean', 'p_value':'count'})
    rna_hybrid_crna_data.rename(columns={'p_value':'off_targets_count_transcripts'}, inplace=True)

    return rna_hybrid_crna_data
"""Functions for main scripts"""

import pandas as pd
import os
import json

from utilities import ensembl_data, ncbi_data, rnacofold_results, rnafold_mirna_results

from utilities_polycistron import mirna_polycistron_data, mirna_sirna_pairs

"""Load config"""

with open(os.path.join("config.json"), "r") as f:
    config = json.load(f)

compl_dict = config["compl_dict"]


def test_ensembl_data():
    """Check information from ensemble database"""

    gene_name = "PDCD1"
    start_gene, end_gene = ensembl_data(gene_name)

    assert (start_gene == 241849684) & (
        end_gene == 241859094
    ), "Ensemble data isn't correct"


def test_ncbi_data():
    """Check information from ncbi database"""

    ncbi_name = "NM_005018"
    _, _, cds_start, cds_end, _ = ncbi_data(ncbi_name)

    assert (cds_start == 57) & (cds_end == 923), "NCBI data isn't correct"


def test_rnacofold_results():
    """Check results from Vienna RNAcofold"""
    vienna_output_directory = os.path.join("test_data", "vienna")
    rnacofold_input = os.path.join("test_data", "data_for_rnacofold.fasta")
    rnacofold_data = rnacofold_results(rnacofold_input, vienna_output_directory)

    assert (
        rnacofold_data.loc[0]["delta_G_cofold"] == "-33.53"
    ), "RNAcofold isn't correct"


def test_rnafold_mirna_results():
    """Check results from Vienna RNAfold"""
    vienna_output_directory = os.path.join("test_data", "vienna")
    rnafold_input = os.path.join("test_data", "data_for_rnafold.fasta")
    mirnafold_data = rnafold_mirna_results(rnafold_input, vienna_output_directory)

    assert mirnafold_data.loc[0]["delta_G_fold"] == -0.1, "RNAfold isn't correct"


def test_mirna_polycistron_data():
    """Check results for natural polycistron data"""

    file_name = os.path.join("databases", "mirna_polycistrons.fa")
    cluster_name = "MIR23AHG"
    members, _, _, _, _ = mirna_polycistron_data(
        file_name, cluster_name, flank_size=200
    )

    assert members == [
        "MI0000081",
        "MI0000085",
        "MI0000079",
    ], "mirna polycistron data isn't correct"


def test_mirna_sirna_pairs():
    """Check results for sirna pairs"""
    mirna_in_polycistrons_path = os.path.join(
        "test_data", "mirna_in_polycistrons_test.csv"
    )
    mirna_in_polycistrons = pd.read_csv(mirna_in_polycistrons_path)
    scaffold_path = os.path.join("test_data", "scaffold.txt")
    with open(scaffold_path) as f:
        scaffold = f.read()
    all_sequences_path = os.path.join("test_data", "all_sequences_test.csv")
    all_sequences_csv = pd.read_csv(all_sequences_path)
    all_sequences = []
    for ind in all_sequences_csv.index:
        all_sequences.append(
            [all_sequences_csv.loc[ind, "Name"], all_sequences_csv.loc[ind, "Sequence"]]
        )

    pairs, _ = mirna_sirna_pairs(
        mirna_in_polycistrons, scaffold, all_sequences, type="similar", k=0.4
    )

    assert (
        pairs.loc[0, "mirna_names"] == "hsa-mir-17"
    ), "mirna/sirna pairs isn't correct"

"""Functions for main scripts"""

import numpy as np
import pandas as pd
import os
import json
from Bio import pairwise2
from pyensembl import EnsemblRelease
from Bio import Entrez
from Bio import SeqIO
import subprocess

"""Load config"""

with open(os.path.join("config.json"), "r") as f:
    config = json.load(f)

compl_dict = config["compl_dict"]


def ensembl_data(gene_name):
    """
    Get data about gene position in genome
    from Ensemble databases
    """
    data = EnsemblRelease(109)
    gene = data.genes_by_name(gene_name)[0]
    start_gene = gene.start - 200
    end_gene = gene.end + 200
    return start_gene, end_gene


def ncbi_data(ncbi_name):
    """
    Get data from NCBI databases about transcript
    """
    Entrez.email = "Your.Name.Here@example.org"
    handle = Entrez.efetch(db="nucleotide", id=ncbi_name, rettype="gb", retmode="text")
    seq_record = [seq_record for seq_record in SeqIO.parse(handle, "gb")][0]
    refseq_sequence = str(seq_record.seq)

    features_list = []
    exons = []
    exon_name = 1
    for feature in seq_record.features:
        features_list.append(
            [
                feature.type,
                int(feature.location.start + 1),
                int(feature.location.end),
                feature.qualifiers,
            ]
        )
        if feature.type == "CDS":
            cds_start = int(feature.location.start + 1)
            cds_end = int(feature.location.end)
        if feature.type == "exon":
            exons.append(
                [
                    "exon" + str(exon_name),
                    int(feature.location.start + 1),
                    int(feature.location.end),
                ]
            )
            exon_name += 1

    return refseq_sequence, features_list, cds_start, cds_end, exons


def correct_sequences(mirna_data, refseq_sequence, compl_dict):
    """
    Check that in table correct strand and fix if not
    (miRNA should be complimentary to transcript)
    """
    mirna_data["start_mirna"] = 0
    mirna_data["end_mirna"] = 0
    for ind in mirna_data.index:
        mirna_i = mirna_data.loc[ind]["Sequence"]
        mirna_i = mirna_i.upper().replace("U", "T")
        mirna_r_i = "".join([compl_dict[s] for s in mirna_i][::-1])
        if mirna_i[:-2] in refseq_sequence:
            mirna_correct = mirna_r_i
            mirna_rev = mirna_i
        else:
            mirna_correct = mirna_i
            mirna_rev = mirna_r_i

        mirna_data.loc[ind, "Sequence"] = mirna_correct

        start_pos = len(refseq_sequence.split(mirna_rev)[0])
        end_pos = start_pos + len(mirna_i)
        mirna_data.loc[ind, "start_mirna"] = start_pos
        mirna_data.loc[ind, "end_mirna"] = end_pos

    mirna_data = mirna_data[pd.notna(mirna_data["Sequence"])]
    return mirna_data


def data_for_vienna(mirna_data, refseq_sequence, file_name, ncbi_name):
    """
    Prepare all neccesary data for Vienna RNA tool
    """
    all_alignments = []
    for _, sirna_i in enumerate(mirna_data["Sequence"]):
        sirna_i = sirna_i.upper().replace("U", "T")
        sirna_r_i = "".join([compl_dict[s] for s in sirna_i][::-1])
        if sirna_i[:-2] in refseq_sequence:
            all_aligh = pairwise2.align.localms(refseq_sequence, sirna_i, 2, -3, -5, -2)

            sirna_start = np.argmax([l != "-" for l in all_aligh[0].seqB])
            all_alignments.append(
                [
                    sirna_i,
                    sirna_r_i,
                    refseq_sequence[sirna_start : sirna_start + len(sirna_i)],
                ]
            )

        elif sirna_r_i[2:] in refseq_sequence:
            all_aligh_r = pairwise2.align.localms(
                refseq_sequence, sirna_r_i, 2, -3, -5, -2
            )

            sirna_start = np.argmax([l != "-" for l in all_aligh_r[0].seqB])
            all_alignments.append(
                [
                    sirna_i,
                    sirna_i,
                    refseq_sequence[sirna_start : sirna_start + len(sirna_i)],
                ]
            )

        else:
            all_aligh = pairwise2.align.localms(refseq_sequence, sirna_i, 2, -3, -5, -2)
            all_aligh_r = pairwise2.align.localms(
                refseq_sequence, sirna_r_i, 2, -3, -5, -2
            )
            if all_aligh[0].score >= all_aligh_r[0].score:
                sirna_start = np.argmax([l != "-" for l in all_aligh[0].seqB])
                all_alignments.append(
                    [
                        sirna_i,
                        sirna_r_i,
                        refseq_sequence[sirna_start : sirna_start + len(sirna_i)],
                    ]
                )
            else:
                sirna_start = np.argmax([l != "-" for l in all_aligh_r[0].seqB])
                all_alignments.append(
                    [
                        sirna_i,
                        sirna_i,
                        refseq_sequence[sirna_start : sirna_start + len(sirna_i)],
                    ]
                )

    all_alignments_results = pd.DataFrame(
        all_alignments, columns=("sirna", "sirna_s", "gene_seq")
    )
    alignments_path = os.path.join(
        "files", "outputs", "vienna", f"{file_name}_for_rnacofold.csv"
    )
    all_alignments_results.to_csv(alignments_path, index=False)

    sequences_path = os.path.join(
        "files", "outputs", "vienna", f"{file_name}_for_rnacofold.fasta"
    )
    with open(sequences_path, "w") as f:
        for _, pair in enumerate(all_alignments):
            f.write("> " + pair[0] + "\n")
            f.write(pair[1] + "&" + pair[2] + "\n")

    transcript_path = os.path.join(
        "files", "outputs", "vienna", f"{file_name}_transcript.fasta"
    )
    with open(transcript_path, "w") as f:
        f.write("> " + ncbi_name + "\n")
        f.write(refseq_sequence + "\n")

    mirna_seq_path = os.path.join(
        "files", "outputs", "vienna", f"{file_name}_for_rnafold_mirna.fasta"
    )
    with open(mirna_seq_path, "w") as f:
        for step, pair in enumerate(all_alignments):
            f.write("> " + pair[0] + "\n")
            f.write(pair[1] + "\n")

    return sequences_path, alignments_path, transcript_path, mirna_seq_path


def rnacofold_results(rnacofold_input, vienna_output_directory):
    """Launch RNAcofold tool"""
    rnacofold_output = rnacofold_input.replace(".fasta", ".txt")
    with open(rnacofold_input, "r") as infile, open(rnacofold_output, "w") as outfile:
        subprocess.run(
            ["RNAcofold", "-p"],
            stdin=infile,
            stdout=outfile,
            check=True,
            cwd=vienna_output_directory,
        )

    with open(rnacofold_output) as f:
        rnacofoldresults = f.readlines()

    rnacofold_data = []
    for i in range(0, len(rnacofoldresults), 5):
        rnacofold_data.append(
            [
                rnacofoldresults[i].split(" ")[1].split("\n")[0],
                rnacofoldresults[i + 4].split("delta G binding=")[1].split("\n")[0],
            ]
        )

    rnacofold_data = pd.DataFrame(
        rnacofold_data, columns=("Sequence", "delta_G_cofold")
    )
    rnacofold_data.drop_duplicates(subset="Sequence", inplace=True)

    return rnacofold_data


def rnafold_results(
    mirna_data,
    all_alignments_results,
    transcript_path,
    vienna_output_directory,
    ncbi_name,
):
    """Launch RNAfold"""
    rnafold_output = transcript_path.replace(".fasta", ".out")
    with open(transcript_path, "r") as infile, open(rnafold_output, "w") as outfile:
        subprocess.run(
            ["RNAfold", "-p", "-d2"],
            stdin=infile,
            stdout=outfile,
            check=True,
            cwd=vienna_output_directory,
        )

    with open(os.path.join(vienna_output_directory, f"{ncbi_name}_dp.ps")) as f:
        rnafold_mrna = "".join(f.readlines())

    mrna_seq = (
        rnafold_mrna.split("sequence { (\\\n")[1]
        .split("\\\n) }")[0]
        .replace("\\\n", "")
        .replace("U", "T")
    )
    bpp = rnafold_mrna.split("%start of base pair probability data\n")[1].split(
        "\nshowpage"
    )[0]
    bbp = pd.DataFrame(
        [i.split(" ") for i in bpp.split("\n")], columns=("p1", "p2", "prob", "type")
    )

    bbp = bbp[bbp["type"] == "ubox"]
    bbp["prob"] = bbp["prob"].astype("float")

    bbp["p1"] = bbp["p1"].astype("int") - 1
    bbp["p2"] = bbp["p2"].astype("int") - 1

    rnafold_data = []
    for sirna_i in mirna_data["Sequence"]:
        gene_s = all_alignments_results[all_alignments_results["sirna"] == sirna_i][
            "gene_seq"
        ].max()
        start_sirna = len(mrna_seq.split(gene_s)[0])
        end_sirna = start_sirna + len(sirna_i)
        p1 = bbp[((bbp["p1"] >= start_sirna) & (bbp["p1"] <= end_sirna))]
        p2 = bbp[((bbp["p2"] >= start_sirna) & (bbp["p2"] <= end_sirna))]
        pair_prob = pd.concat([p1, p2])["prob"].mean()
        rnafold_data.append([sirna_i, pair_prob])

    rnafold_data = pd.DataFrame(rnafold_data, columns=("Sequence", "pair_probability"))
    rnafold_data.drop_duplicates(subset="Sequence", inplace=True)
    return rnafold_data


def rnafold_mirna_results(mirna_seq_path, vienna_output_directory):
    """Launch RNAcofold tool for mirna"""
    rnafold_output = mirna_seq_path.replace(".fasta", ".out")
    with open(mirna_seq_path, "r") as infile, open(rnafold_output, "w") as outfile:
        subprocess.run(
            ["RNAfold", "-p", "-d2"],
            stdin=infile,
            stdout=outfile,
            check=True,
            cwd=vienna_output_directory,
        )

    with open(rnafold_output) as f:
        mirnafold = f.readlines()

    mirnafold_data = []
    for i in range(0, len(mirnafold), 6):
        mirnafold_data.append(
            [
                mirnafold[i].split(" ")[1].replace("U", "T").split("\n")[0],
                float(mirnafold[i + 2].split(" ")[-1].split(")")[0].replace("(", "")),
            ]
        )

    mirnafold_data = pd.DataFrame(mirnafold_data, columns=("Sequence", "delta_G_fold"))
    mirnafold_data.drop_duplicates(subset="Sequence", inplace=True)

    return mirnafold_data


def mrna_region(start, end, cds_start, cds_end):
    """Get transcript region"""
    if (start >= cds_start) & (end <= cds_end):
        return "CDS"
    elif start < cds_start:
        return "5-UTR"
    else:
        return "3-UTR"


def exon_name(start, exons):
    """Get exon name"""
    return list(filter(lambda p: (start >= p[1]) & (start <= p[2]), exons))[0][0]


def transcript_features(mirna_data, cds_start, cds_end, exons):
    """Collect information about mirna/sirna"""
    mirna_data["region"] = 0
    mirna_data["exon"] = 0
    for ind in mirna_data.index:
        start_mirna = mirna_data.loc[ind]["start_mirna"]
        end_mirna = mirna_data.loc[ind]["end_mirna"]
        region = mrna_region(start_mirna, end_mirna, cds_start, cds_end)
        exon = exon_name(start_mirna, exons)
        mirna_data.loc[ind, "region"] = region
        mirna_data.loc[ind, "exon"] = exon
    return mirna_data


def compute_GC_context(sequence):
    """GC context for mirna/sirna"""
    return (list(sequence).count("G") + list(sequence).count("C")) / len(sequence)

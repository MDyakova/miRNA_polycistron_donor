"""
Functions to prepare Gene Bank file for SnapGene tool
"""

import os

"""
Next templates for gene bank file.
This file allow to automatically make SnapGene file with all features and primers
and open it in SnapGene Viewer without pair SnapGene program.
"""

TITLE = """LOCUS       {gene_name}        {len_full_sequence} bp DNA     linear   UNA {date_today}
DEFINITION  {gene_name}.
ACCESSION   .
VERSION     .
KEYWORDS    .
SOURCE      synthetic DNA construct
  ORGANISM  synthetic DNA construct
REFERENCE   1  (bases 1 to {len_full_sequence})
  AUTHORS   .
  TITLE     Direct Submission
  JOURNAL   For SnapGene Viewer
            https://www.snapgene.com
FEATURES             Location/Qualifiers"""

FEATURE_SOURCE = '''     source          1..{len_full_sequence}
                     /mol_type="other DNA"
                     /note="color: #ffffff"
                     /organism="synthetic DNA construct"'''

ORIGIN = """ORIGIN
{origin_seq}
//"""


def misc_feature_template(start, end, label, color, direction):
    """
    Make a template for snapgene features
    """

    if direction == "None":
        misc_f = f'''     misc_feature    {start}..{end}
                        /label={label}
                        /note="color: {color}"'''
    else:
        misc_f = f'''     misc_feature    {start}..{end}
                        /label={label}
                        /note="color: {color}; direction: {direction}"'''
    return misc_f


def primer_template(name, seq, date_today, start, end):
    """
    Make a template for snapgene primers
    """
    primer_f = f'''     primer_bind     complement({start}..{end})
                     /label={name}
                     /note="color: black; sequence:
                     {seq}; added:
                     {date_today}"'''
    return primer_f


def gene_bank_file(
    gene_name,
    full_sequence,
    date_today,
    elements_list,
    files_name,
    output_folder,
    oligos=None,
    title=TITLE,
    feature_sourse=FEATURE_SOURCE,
    origin=ORIGIN,
):
    """
    Make gene bank file for snapgene tool
    """

    if oligos is None:
        oligos = []

    title = title.format(
        gene_name=gene_name,
        full_sequence=full_sequence,
        date_today=date_today,
        len_full_sequence=len(full_sequence),
    )

    feature_sourse = feature_sourse.format(len_full_sequence=len(full_sequence))

    all_misc_feature = ""
    for feature in elements_list:
        start = feature[1]
        end = feature[2]
        name = feature[0]
        color = feature[4]
        direction = feature[3]
        if direction == "+":
            direction = "RIGHT"
        elif direction == "-":
            direction = "LEFT"
        else:
            direction = "None"
        misc_feature = misc_feature_template(start, end, name, color, direction)
        all_misc_feature += misc_feature + "\n"

    all_primers = ""
    for oligo in oligos:
        start = oligo[2]
        end = oligo[3]
        name = oligo[0]
        seq = oligo[1]

        primer_feature = primer_template(name, seq, date_today, start, end)
        all_primers += primer_feature + "\n"

    # This part make a sequence in Gene Bank format.
    origin_seq = ""
    for i in range(len(full_sequence)):
        if i % 60 == 0:
            origin_seq += "\n" + " " * (9 - len(str(i + 1))) + str(i + 1) + " "
        elif i % 10 == 0:
            origin_seq += " "
        origin_seq += full_sequence[i].lower()
    origin_seq = origin_seq[1:]

    origin = origin.format(origin_seq=origin_seq)

    gbk_file_name = os.path.join(output_folder, files_name + ".gbk")

    with open(gbk_file_name, "w", encoding="utf-8") as file:
        file.write(title + "\n")
        file.write(feature_sourse + "\n")
        file.write(all_misc_feature)
        file.write(all_primers)
        file.write(origin + "\n")


def find_elements(cds_start, cds_end, exons, mirna_data):
    """
    Make lists with CDS and exons positions for SnapGene file.
    """
    elements_list = []

    # Add CDS
    elements_list.append(["CDS", cds_start, cds_end, "None", "#40139c"])

    # Add exons
    for exon in exons:
        name_exon = exon[0]
        start_exon = exon[1]
        end_exon = exon[2]
        elements_list.append([name_exon, start_exon, end_exon, "None", "#e3d914"])

    # Add mirna/sirna
    oligos = []
    for ind in mirna_data.index:
        start_mirna = mirna_data.loc[ind]["start_mirna"]
        mirna_seq = mirna_data.loc[ind]["Sequence"]
        name_mirna = "_".join([mirna_data.loc[ind]["Name"], str(start_mirna)])
        oligos.append([name_mirna, mirna_seq, 1, len(mirna_seq)])

    return elements_list, oligos


def find_elements_polycistron(
    new_scaffold_seq, all_sequences, sirna_names, all_new_sequences
):
    """
    Make lists with sirna/mirna and hairpin sequences for SnapGene file.
    """
    elements_list = []

    for sirna_name, sirna_seq in zip(sirna_names, all_new_sequences):
        sirna_start = len(new_scaffold_seq.split(sirna_seq)[0])
        sirna_end = sirna_start + len(sirna_seq)
        elements_list.append(
            [f"{sirna_name}_hairpin", sirna_start, sirna_end, "None", "#c96f0e"]
        )

    oligos = []
    for sirna_name, sirna_seq in all_sequences:
        if sirna_name in sirna_names:
            sirna_start = len(new_scaffold_seq.split(sirna_seq)[0])
            sirna_end = sirna_start + len(sirna_seq)
            oligos.append([sirna_name, sirna_seq, 1, len(sirna_seq)])
    return elements_list, oligos

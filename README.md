
# miRNA Polycistron Donor Design Tool

This tool facilitates the design of donor sequences containing artificial miRNA polycistrons. It integrates polycistrons into DNA using CRISPR-Cas9 and HDR techniques, enabling the expression of selected miRNA. Additionally, the tool computes critical characteristics of miRNA/siRNA that influence silencing efficiency, providing users with insights to optimize gene regulation and therapeutic applications.

## Why This Tool Is Important

miRNA and siRNA play a pivotal role in post-transcriptional gene regulation. Designing artificial polycistrons allows for precise expression of multiple miRNA/siRNA sequences, which can:
- Enhance gene silencing efficiency in experimental and therapeutic settings.
- Provide a platform to study natural polycistron configurations and compare their functionality.
- Automate complex design processes, making advanced miRNA/siRNA experiments more accessible.

---

## Features

- **Docker-based environment** for easy setup and usage.
- **Automatic computation** of efficiency features for miRNA/siRNA sequences.
- **Integrated folding analysis** to validate polycistron design.
- **Compatibility** with natural polycistrons from the human genome.
- Support for **customized polycistron** scaffolds and default settings.
- Output in **CSV, SnapGene files, and graphical representations**.

---

## Installation Instructions

1. **Clone the Repository:**
   Clone this repository to your local machine:
   ```bash
   git clone https://github.com/MDyakova/miRNA_polycistron_donor.git
   cd miRNA_polycistron_donor
   ```

2. **Set Up the Output Directory:**
   Ensure that the folder structure for outputs and test results is in place:
   ```
   miRNA_polycistron_design_launch/files/data
   miRNA_polycistron_design_launch/files/outputs
   ```
   Save your input files (e.g., miRNA/siRNA table) inside `files/data`.

3. **Pull Docker Container:**
   ```bash
   docker pull mdyakova/mirna_polycistron_donor_design:latest
   ```

3. **Run the Docker Container:**
   Launch the Docker container for the two main steps:

   **Step 1: miRNA Analysis**
   Provide information about your mirna/sirna.
   ```bash
   docker run -v ${PWD}/miRNA_polycistron_design_launch/files:/main/files -d mirna_polycistron_donor_design:latest /main/miRNA_analysis.py
   ```

   **Step 2: Polycistron Construction**
   Make polycistron construct with your mirna/sirna
   ```bash
   docker run -v ${PWD}/miRNA_polycistron_design_launch/files:/main/files -d mirna_polycistron_donor_design:latest /main/polycistron_construct.py
   ```
---

## Usage Guidelines

### Step 1: Prepare Input Data
1. Save your table of miRNA/siRNA sequences in the `miRNA_polycistron_design_launch/files/data` directory.
2. The table should have the following format:
   | Name      | Sequence               | Source        |
   |-----------|------------------------|---------------|
   | miRNA1    | AAGGTTCCGGAACTCCATCC   | Literature    |
   | miRNA2    | GGCCTTAACGGTTCCAAGGG   | Database      |
   | miRNA3    | TTCGGAAACCGGTTCCAATC   | Experimental  |

   Include at least six miRNA/siRNA sequences for natural polycistron with five members.
   If you do use the default scaffold, please ensure that the number of sequences is greater than the number of members + 1.

### Step 2: Update the Configuration File
Edit the `input_data.json` file inside the "miRNA_polycistron_design_launch" folder. Update the fields with your project-specific details:
```json
"general_information": {
    "file_name" : "{your_file_name}",
    "output_folder_name":"{your_output_file_name}"
},
"gene_information": {
    "gene_name" : "{gene name}",
    "ncbi_name" : "{transcript id from ncbi}",
    "chromosome_name" : "{chromosome name for your gene}"
},
"data" : {
    "mirna_table" : "files/data/{file with your mirna/sirna}"
}
```

### Step 3: Review miRNA/siRNA Efficiency Features

After running the miRNA analysis script, open the output file with the suffix `_mirna_with_features.csv`. This file contains important features to help you select the best miRNA/siRNA. Users should choose the miRNA/siRNA by setting the `choice` column to `1` for the selected sequences. Here is an example of what the table might look like:

| Name       | Sequence              | Region  | Exon Name | GC Content | ΔG Target Binding | ΔG Folding | Pair Probability | Choice |
|------------|-----------------------|---------|-----------|------------|-------------------|------------|------------------|--------|
| miRNA1     | AAGGTTCCGGAACTCCATCC  | 3'-UTR  | Exon_3    | 45%        | -20.5 kcal/mol    | -5.2 kcal/mol | 0.85            | 0      |
| miRNA2     | GGCCTTAACGGTTCCAAGGG  | CDS     | Exon_5    | 50%        | -18.7 kcal/mol    | -4.8 kcal/mol | 0.78            | 1      |
| miRNA3     | TTCGGAAACCGGTTCCAATC  | 5'-UTR  | Exon_1    | 55%        | -22.3 kcal/mol    | -6.0 kcal/mol | 0.90            | 0      |

### Explanation of Columns:
- **Name:** Identifier for the miRNA/siRNA.
- **Sequence:** The nucleotide sequence.
- **Region:** The region of the transcript the sequence targets (3'-UTR, CDS, or 5'-UTR).
- **Exon Name:** Name or index of the exon containing the target site.
- **GC Content:** Percentage of G and C bases in the sequence, which influences stability.
- **ΔG Target Binding:** Free energy of binding to the target transcript (lower values indicate stronger binding).
- **ΔG Folding:** Free energy of miRNA/siRNA folding (lower values indicate potential self-complementary structures to avoid).
- **Pair Probability:** Probability of transcript folding in the target region (higher values indicate less accessibility).
- **Choice:** Set to `1` for selected miRNA/siRNA and `0` for others.

---

This table will help users make an informed decision by selecting miRNA/siRNA with better silencing efficiency.

### Step 4: Run Polycistron Construction
Once you’ve marked your selected miRNA/siRNA, execute the polycistron construction script and check the results in:
```
files/data/outputs/{output_folder_name}
```
---

## Results
All output files will be located in:
```
files/data/outputs/{output_folder_name}
```
This includes:
- Processed miRNA/siRNA data.
- Polycistron designs.
- RNA folding analysis.
- CSV and SnapGene files for further use.

---

## Requirements

The following dependencies are installed automatically within the Docker container:

- `pandas==1.5.3`
- `pyensembl==2.2.8`
- `matplotlib==3.7.1`
- `biopython==1.79`
- `seaborn==0.12.2`
- `openpyxl==3.1.2`
- `pytest==6.2.5`
- `ViennaRNA-2.7.0`

> **Note:** The ViennaRNA package is specifically used for RNA folding computations.

---

## License

This tool incorporates **ViennaRNA** for RNA folding computations. For compliance, mention ViennaRNA in your documentation and dependency file if applicable:
> This tool uses ViennaRNA for RNA folding calculations. For more information, visit [ViennaRNA](https://www.tbi.univie.ac.at/RNA/).

---

Let me know if additional edits are needed!
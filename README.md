# MODesign – iSBH‑sgRNA Design Toolkit

MODesign is a four‑script pipeline that helps you design modular iSBH‑sgRNA sensors, rank the best candidates and generate cloning‑ready iSBH-sgRNA oligos. These scripts are associated with the publication:

# Specific Modulation of CRISPR Transcriptional Activators through RNA-Sensing Guide RNAs in Mammalian Cells and Zebrafish Embryos
 Oana Pelea, Sarah Mayes, Quentin R.V. Ferry, Tudor A. Fulga, Tatjana Sauka-Spengler
eLife. https://doi.org/10.7554/eLife.87722.2 

## Scripts

**Script_1a_MODesign_designing_spacer_star_sequences.py**  
Finds spacer / spacer* pairs that fold into the desired mini‑hairpin. Run first.

**Script_1b_MODesign.py**  
Combines a chosen spacer pair with your RNA trigger(s) to build modular iSBH‑sgRNA designs. Code also gives you the option to select desired sizes for the iSBH-sgRNA loop.

**Script_2_MODesign_scoring_algorithm.py**  
Scores and ranks designs so you can pick the most promising ones for wet‑lab testing.

**Script_3_automating_ISBH_oligo_design.py**  
Outputs cloning oligos compatible with Addgene vectors 200234 (iSBH-sgRNA expression from U6 promoters), 200235 (small triggers complementary with loop and spacer * sequences to be expressed from a mammalian U6 promoter) & 200237 (1xCTS sequences for fluorescent reporters) and simulates final expected plasmid maps (after cloning).

Example input files are provided for each script; please use the exact same format and just change RNA sequences and sequence names with your custom inputs.

## Prerequisites

| Requirement | Version | Notes |
| --- | --- | --- |
| Python | 2.7.x | Code won’t run under Python ≥3.0 without edits. |
| NUPACK | 3.2.2 | Required for RNA folding; make sure the nupack binaries are in your $PATH. Please download this from [https://www.nupack.org](https://www.nupack.org). |

## Contact

For questions or feedback, email Oana Pelea: o.pelea@bioc.uzh.ch or oanapelea1@gmail.com

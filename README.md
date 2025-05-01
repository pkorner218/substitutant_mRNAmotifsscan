# substitutant_mRNAmotifsscan
Take a library of substitant peptides from mass spectronomy and scan for mRNA motifs around the substitution site


## Summary
This small pipeline was created to provide a toolset that allows the analysis of nucleotide and amino acid frequencies surrounding substitutant sites detected in mass spectronomy data.

## Input

A list or file '(first column is ID) of proteomics IDs of Substitutant peptides following the scheme:
"ENSTID_Genename_Codon_Substitutantposition-within-peptide_Substututant(H2Q)_Peptidestart_Peptideend_Peptidesequence"
"ENST00000XXXXXX_XXX_CAC_9_H2Q_1_15_AAAGAAATQLEVAR" 

A list or file '(first column is ID) of proteomics IDs of WT peptides following the scheme:
"ENSTID_Genename_Codon_WT_Peptidestart_Peptideend_Peptidesequence"
"ENST00000XXXXXX_XXX_CAC_WT_549_568_AAISSGIEDPVPTLHLTER"

## Technical

The pipeline is a set of scripts in python, partially utilizing R for plotting.
In the following several of the main packages which are required are listed.
* python3
* R
* pandas
* numpy
* pheatmap

## Usage

`python3 calculate_signature_frequencies.py -i [infile] -o [outfoldername] -s [substitution] -r [reference file] -t [type]`

This script takes the list, maps all peptides (either WT or given Substitutant) to a given reference file and creates all frequency files in the outfolder. It works on substitutant peptides or WT peptides based on the -t flag (S or W). e.g.

`python3 calculate_signature_frequencies.py -i eH2Q.txt -o "./Sub" -s HQ -r ../New_CDS_Framevalidated_gencode.v44.fa -t S`\
`python3 calculate_signature_frequencies.py -i eWT.txt -o "./WT" -s HQ -r ../New_CDS_Framevalidated_gencode.v44.fa -t WT`

The following script takes two of the created frequencyfiles and plots their difference as a lineplot.

`python3 plot_signature_lines.py -is [substitutant infile] -iwt [WT infile] -o [name of outfile pdf]`\
`python3 plot_signature_lines.py -is ./Sub/HQ_Anuc.txt -iwt ./WT/HQ_Anuc.txt -o ./lineplots/A_nuc_difference.pdf`

The following script automatically creates nucleotide and codon frequency heatmap plots based on given directories.

`python3 calculate_signature_heatmap.py -i1 [folder with substitutant freuqncy files] -i2 [folder with WT frequency files] -o [output folder]`\
`python3 calculate_signature_heatmap.py -i1 ./Sub -i2 ./WT -o ./heatmap_plots/`

## heatmap plots in R

Heatmap plots were created using R. The automated python results of "calculate_signature_heatmap.py" can be read in R and plotted by the following example lines. 

`df <- read.table("CAC_Codon_frequency_difference.txt", sep = "\t", header = T, row.names = "nucleotide", check.names=FALSE)`\
`pheatmap(df, cluster_col = F, cluster_row = F, border_color = "black", cellheight = 30, cellwidth = 30, breaks = seq(-0.2,0.2, length.out=(100)))` 








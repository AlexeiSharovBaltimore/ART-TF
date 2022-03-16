ART-TF project includes programs for bioinformatic analysis of ChIP-seq data sets and comparison with 
gene expression data.

1. combine_peaks_TFs.pl

Syntax: combine_peaks_TFs.pl file_list.txt peaks_all.txt -TSS TSS_symbols_hg19.txt

Perl script to associate ChIP-seq peaks with genes and assemble data from many ChIP-seq experiments into a single file.
File file_list.txt is input tab-delimited text file with first column containing names of ChIP-seq experiments, and
second column containing file names of corresponding ChIP-seq results. Each ChIP-seq result file is a tab-delimited text file
with 4 columns: chromosome, peak start position, peak end position, peak score. Chromosome is formatted as "chr1". File
TSS_symbols_hg19.txt is a tab-delimited text file with 8 columns: gene symbol, chromosome (e.g., "chr1"), position of TSS,
strand (0 = forward, 1 = reverse), length of transcript (mRNA), number of exons, and accession number/ID.

2. generate_targets.pl

Syntax: generate_targets.pl peaks_all.txt TF_targets_human_add_far.txt -o 1 -add

Perl script to generate gene sets containing target genes for various transcription factors (TFs). 
Input file peaks_all.txt is the output of program
combine_peaks_TFs.pl. Output file TF_targets_human_add_far.txt is a tab-delimited text file formatted as a gene-set file
for ExAtlas. Each gene set is described by a line where the first item is the name of a gene set, the second item is a 
description or other characteristic, such as gene count, then go gene symbols (from the third item to the end of the 
line). Attributes of genes are optional and require additional lines. In particular, here the output file includes score
as the attribute. Attribute lines in the file start wth a tab, which means that the first item is empty, and the second 
item is the name of attibute (here it is "score"). Scores for each gene in the set are specified starting from the third
item and to the end of the line. The number of score values should match the number of genes in the previous line.

Options for generate_targets.pl are: option "-o 1" means that only distal binding cites are analyzed (at distances 
from 0.5 Kb to 100 Kb from TSS on both sides); option "-o 2" means that only proximal binding cites are analyzed (at distances 
from -0.5 Kb to +0.5 Kb from TSS); if option "-o" is not used, then all binding sites are analyzed. Option "-add" means
that ChIP-seq scores are added for all binding sites linked to each gene. Each ChIP-seq peak was associated with a 
maximum of 3 genes whose TSS was within 100 Kb from the peak center. If th program is run without option "-add", then 
only one binding site is considered with the maximum score. Scores of gene/peak associations were calculated 
as symbol quality multiplied by the binding score (ChIP-seq) and divided by the distance from the peak to TSS 
(Kb, capped at 1 Kb). Symbol quality was equal 1 for "weak" symbols (e.g., containing 4 digits in a row, or strings 
"FAM", "MIR", "MRP", and "orf") and 3 for normal symbols. Genes with association scores <20% of the maximum value 
(i.e., for the best matching gene) were not reported as associated with the given ChIP-seq peak. 

3. The next step of the analysis is finding regulated targets of TFs; it is done using ExAtlas after gene 
sets of target genes are uploaded. You don't need to
upload them again because these files are available in ExAtlas as public resources with the following names:
public-ART-TF_taergets_hg19_enhancers_add and public-ART-TF_taergets_hg19_promoters_add, which specify traget genes
based on distal binding sites and proximal binding sites, respectively. ExAtlas is a software for meta-analysis
of gene expression data. Besides standard statistical analysis of gene expression (similar to NIA Array Analysis) 
it supports several methods for meta-analysis and generates results for all combinations of data in 
multi-component data sets (e.g., all gene expression profiles vs. all GO annotations). ExAtlas demo is available 
at http://alexei.nfshost.com/exatlas/ (Full installation of ExAtlas: http://webtools.systemsmedicine.jp/exatlas/) and the code can be accessed at 
https://github.com/AlexeiSharovBaltimore/ExAtlas. To enter as guest, click button "Start using ExAtlas". Then, to access 
ART-TF data, users need to select organism species by clicking the pull-down menu "Select organism" shown by blue 
arrow at the top of the screen. Please, select "Human (Homo sapiens)".

To identify regulated targets of TFs, locate section "Select Data Files" and click on the pull-down menu in the first
line (gene expression profiles). Select data "public-CREST_Human_TF_induction_ESC_20181203" which is large-scale
study of the gene expression change in ES cells after induction of 510 individual TFs (Nakatake et al. 2020. Cell Reports,
31(7):107655) and click button "Open" in the same line at the left. In the next screen, locate section 4 "Other functions"
and click on "Geneset analysis ...". In the next screen select geneset file "public-ART-TF_taergets_hg19_enhancers_add"
and check the box "Use gene attributes". This checkbox is needed to instruct the program to use gene attributes, which
in our case is the score of TF binding to enhancer(s) of each target gene. Also check another checkbox "Identify 
associated genes", and then click the button "Geneset analysis". The analysis may take 25-30 minutes; click button
"Check your task" when it is available; also you can click in "Log file" to see the progression.

4. parse_targets.pl

Syntax: parse_targets.pl TF_regulated_targets_EPFP03_add_far.txt TF_regulated_targets_EPFP03_add_near.txt target_scores_add.txt -ind indirect_list_new.txt

Perl script to parse output files with regulated gene targets (see step #3) downloaded from ExAtlas and named 
"TF_regulated_targets_EPFP03_add_far.txt" (for enhancers) and "TF_regulated_targets_EPFP03_add_near.txt" (for promoters).
Option " -ind indirect_list_new.txt" is needed to supply the program with a file that specifies acceptable surrogate TFs
so that ChIP-seq surrogate data can be used in case there is no ChIP-seq data for some induced TF. Surrogate ChIP-seq data
also makes lists of regulated target genes more complete to increase their explanatoiry power.
The output file target_scores_add.txt is a tab-delimited text file with 11 columns: (a) ChIP-seq data name, (b) Induced TF,
(c) direct=1 or indirect=0 type of regulated target genes (indirect means using surrogate ChIP-seq data), (d) z-value for
gene set enrichment of target genes with enhancer binding among upregulated genes, (e) z-value for gene set enrichment 
of target genes with promoter binding among upregulated genes, (f) z-value for gene set enrichment of target genes with 
enhancer binding among downregulated genes, (g) z-value for gene set enrichment of target genes with promoter binding 
among downpregulated genes, (h) number of upregulated target genes with enhancer binding, (i) number of upregulated 
target genes with promoter binding, (j) number of downregulated target genes with enhancer binding, (k) number of 
downregulated target genes with promoter binding.

5. output2targets.pl

Syntax: output2targets.pl TF_regulated_targets_EPFP03_add_far.txt TF_regul_targets_add_far.txt -ind indirect_list_new.txt

Perl script to extract gene sets of regulated gene targets from output files with regulated gene targets (see step #3)
downloaded from ExAtlas and named "TF_regulated_targets_EPFP03_add_far.txt" (for enhancers). Option 
" -ind indirect_list_new.txt" is needed to supply the program with a file that specifies acceptable surrogate TFs.
Output file TF_regul_targets_add_far.txt specifies gene symbols for regulated target genes generated from individual
ChIP-seq experiments and one attribute: EPFP (expected proportion of false positives). 

6. target_overlap.pl

Syntax: target_overlap.pl TF_regul_targets_add_far2021.txt target_overlap_add_far.txt

Syntax: target_overlap.pl TF_regul_targets_add_far2021.txt TF_regul_targets_add_far_combined.txt -g

Perl script to count overlapping regulated gene targets for pairs of gene sets obtained from different ChIP-seq
experiments with the same TF gene. If option "-g" is added at the end, then the program generates combined gene sets
of regulated target genes from multiple ChIP-seq experiments with the same TF gene. For each TF, 4 gene sets can be
generated: upregulated direct targets, downregulated direct targets, upregulated indirect targets, and downregulated
indirect targets. Two attributes are: smallest EPFP among ChIP-seq data sets, and the number of ChIP-seq data sets that
support each regulated target gene. If 3 or more ChIP-seq data sets are available for a TF, then regulated target 
genes suported by a signle ChIUP-seq data set are not included into a combined gene set of regulated targets.

7. compare_indirect_direct.pl

Syntax: compare_indirect_direct.pl target_overlap_add_far2021.txt freq_indirect_direct_far2021.txt

Perl script to compare the quality of direct and indirect regulated targets. The significance of overlap between sets 
of direct and indirect regulated targets for the same TF was quantified by the hyporgeometric test (z-value). The 
overlap between sets of direct and indirect regulated targets was generally lower than the overlap between sets of 
direct regulated targets identified using different ChIP-seq data for the same TF, as follows from the probability 
distribution of z-values.

8. target_far_near_new.pl

Syntax: target_far_near_new.pl TF_regul_targets_add_far_combined2021.txt TF_regul_targets_add_near_combined2021.txt regul_targets_table_add_2021.txt

Syntax: target_far_near_new.pl TF_regul_targets_add_far_combined2021.txt TF_regul_targets_add_near_combined2021.txt regul_all_targets_genesets_add.txt -gall -anova anova-TFs_Sep2018_all.txt -coord TF_targets_human_add_far2018.txt TF_targets_human_add_near2018.txt -indir indirect_list_new.txt

Perl script to generate a table of varyoius categories of regulated gene counts: upregulated and downregulated, direct and indirect,
bound by TFs in enhancers or prtomoters or both, and significantly changed expression (>2 fold and FDR<0.05). Options "-gall",
"-gdir", and "-gind" are used to generate a file with genesets of regulated targtes supported by all ChIP-seq data, non-surrogate
ChIP-seq data, and surrogate ChIP-seq data, respectively. These gene sets are characterized by the following attributes:
(a) smallest EPFP, (b) average logratio (log10) of gene expression channge after induction of the TF, (c) counts of ChIP-seq
data sets supporting the target gene. and (d) genome position of the best binding site of the TF (hg19) or a pair of binding sites
one in enhancer and another in promoter.

9. count_2fold_targets.pl

Syntax: count_2fold_targets.pl regul_all_targets_genesets_add.txt signif_Nov2017_2fold_MedEmerCag_genesets.txt regul_all_targets_2fold.txt

Perl script to generate genesets of regulated targets that also show signiuficant change of expression (>2 fold & FDR<0.05)
in ES cells after induction of corresponding TFs. Both input files and the output file (regul_all_targets_2fold.txt) are
gormatted as genesets for ExAtlas (as explained above).

10. generate_table_targets_v2.pl

Syntax: generate_table_targets_v2.pl regul_dir_targets_genesets_add.txt regul_all_targets_genesets_add.txt output.txt

Perl script to generate a summary table of regulated target genes for all TFs. Output file is a tab-delimited table/text
with 11 columns. Column headers: a) Induced transcription factor (TF) gene symbol; b) Direction of gene expression 
change after TF induction; c) Target gene number ordered by increasing EPFP; d) Target gene symbol; e) Log10-ratio of
target gene expression change after TF induction; f) All regulated targets (including surrogate ChIP-seq): EPFP; 
g) All regulated targets: Number of supporting ChIP-seq data sets; h) All regulated targets: genome position (center); 
i) Direct regulated targets (no surrogate ChIP-seq): EPFP; j) Direct regulated targets: Number of supporting ChIP-seq
data sets; k) Direct regulated targets: genome position (center).




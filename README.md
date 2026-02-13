### 
Marine fungi harbor, and are co-diversified with, remnants of endogenous giant viral elements. 
---
#### By: Syrena Whitner 
---
### To cite this work or code 
#### Whitner, S., Amend, A. S., Culley, A. I. (2026). Marine fungi harbor, and are co-diversified with, remnants of endogenous giant viral elements. [Manuscript submitted for publication to ...].
#### Whitner, S. (2025). Marine fungi harbor, and are co-diversified with, remnants of endogenous giant viral elements. Zenodo. 10.5281/zenodo.18626390.
---
### Abstract
Nucleocytoplasmic large DNA viruses (NCLDVs) are ancient associates of eukaryotes whose endogenized remnants are increasingly recognized as drivers of genome evolution. In this study, we identified NCLDV-like endogenous viral elements in 18 marine Ascomycota fungi isolates from Hawaiian coastal waters. Across viral-bearing contigs, we recovered a diverse set of canonical NCLDV hallmark genes co-localized with auxiliary metabolic and regulatory genes (AMGs). These AMGs, including glycosyltransferases, serine/threonine and F10-like kinases, and short-chain dehydrogenases, suggest that the ancestral viruses encoded pathways influencing carbon metabolism and redox balance in their fungal hosts. The additional recovery of genes such as mRNA capping enzymes and ribonucleotide reductases suggests these viruses were once highly autonomous and contained machinery to support replication independent of the fungal host. The intronization of these genes within fungal scaffolds further supports the potential functional assimilation following endogenization. Phylogenetic and cophylogenetic analyses place these elements as sister to the Mycodnaviridae, supporting a long-standing shared evolutionary history between fungi and NCLDV viruses and pointing to an unrecognized contribution of these viral lineages to fungal genome diversification.

### Script information:
#### Scripts should generally be run in the following order
1. Genomad_CandidateContigs.sh - Initial screening of fungal genomes for viral sequences using GENOMAD (v1.11.2), and then additional contig curation based on values such as gv_marker, taxonomic classification, etc. 
2. Viral_HMM.sh - Predicting and annotating viral / fungal proteins, classification and extraction of viral marker genes
3. Viral_Fungal_Phylogeny.sh - aligns, trims, and constructs 1) fungal host phylogeny and 2) viral phylogeny 
4. parse_augustus_prodigal.sh - parses output files from Augustus (fungal) and prodigal (viral), used for downstream analysis in Contig_plots.R
5. PairwiseSoftAlignments_analysis.sh - Protein embedding based similarity analysis, compile protein sequences from study isolates, NCLDV references, and fungal homologs, generate protein embeddings using large language model–based protein encoders (Harrigan et al. 2024), and calculate pairwise similarity scores for all protein combinations. output used in downstream analysis in PCoA_SoftAlign.R
6. Contig_plot.R - generates contig overlay plots, FIGURE 2
7. PCoA_SoftAlign.R - generates PCoA plots comparing sequences in this study with reference NCLDVs and fungal homologs, FIGURE 3
8. PACo_analysis.R - procrustean approach to cophylogeny analysis, generates tanglegram based on relationship between fungal/viral phylogenies, FIGURE 4
9. Heatmap.R - generates heatmap denoting number of hits for various GVOGs, SUPPLEMENTARY FIGURE 1. 

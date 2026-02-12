########################################################################
#
# Protein embedding based similarity analysis - PART 1
# Compile protein sequences from study isolates, NCLDV references, and fungal homologs
# Generate protein embeddings using large language model–based protein encoders (Harrigan et al. 2024)
# Calculate pairwise similarity scores for all protein combinations
# Convert similarity scores into a distance matrix
# Syrena Whitner 
# December 2025
#
########################################################################

########################################################################
# Step 1:
# pull reference genes from NCBI using esearch 
########################################################################

esearch -db protein -query "DNA topoisomerase II[Title] AND fungi[Organism]" \
  | efetch -format fasta -stop 25 > topoII_fungi_25.faa


genes=(
  "DEAD-box helicase[Title] AND fungi[Organism]"
  "DNA-directed RNA polymerase beta subunit[Title] AND fungi[Organism]"
  "DNA-directed RNA polymerase alpha subunit[Title] AND fungi[Organism]"
  "DNA polymerase family B[Title] AND fungi[Organism]"
  "Transcription initiation factor IIB[Title] AND fungi[Organism]"
  "A32 packaging ATPase[Title] AND viruses[Organism]"
)

names=(SFII RNAPS RNAPL PolB TFIIB A32)

for i in ${!genes[@]}; do
  q="${genes[$i]}"
  n="${names[$i]}"
  echo "Fetching $n ..."
  esearch -db protein -query "$q" \
    | efetch -format fasta -stop 25 > "${n}_25.faa"
done

###
# NOTE : only ended up using the following genes in the analysis: 
# SFII : GVOGm0013
# RNAPS : GVOGm0023
# RNAPL : GVOGm0023 
# TopoII : GVOGm0461
###

########################################################################
# Step 2:
# pull 50 random sequences for each of the above genes from the ICTV list
########################################################################

# i.e just got one gene 
awk '/^>/{f++} f<=50' GVOGm0461_ictv.faa > GVOGm0461_ictv_50.faa

###
# NOTE : at this point, I have 50 random NCLDVs, my sequences, the 25 random fungal references, 
# and the Mycodnoviridae sequences (umich). These are what will be fed into the analysis 

########################################################################
# Step 3:
# concatonate all genes into one files .faa file 
# generate list of taxa names 
########################################################################

# i.e for one gene 
cat *_GVOGm0013.faa > GVOGm0013_cat.faa 

# cleaning up seqs now within all directories just to make everything good to go
find . -type f -name "cat_*.faa" | while read f; do
  echo "Cleaning $f ..."
  awk '/^>/{                                      # Header lines
           gsub(/ .*/, "");                       # 1. remove everything after first space
           gsub(/\|/, "_");                       # 3. replace "|" with "_"
           gsub(/[[:space:]]+$/, "");             # 4. remove trailing spaces
           print; next
       }
       {
           gsub(/\*/, "");                        # 2. remove "*" anywhere in sequence
           gsub(/[[:space:]]+$/, "");             # 4. remove trailing spaces
           print
       }' "$f" > "${f%.faa}_cleaned.faa"
done

# getting all the sequence names 
find . -type f -name "*_cleaned.faa" | while read f; do
  echo "Extracting headers from $f ..."
  out="${f%.faa}_names.txt"
  grep "^>" "$f" | sed 's/^>//' > "$out"
done

### 
# NOTE: at this point I have a list of sequence headers in a .txt file, that I had to manually enter into an excel sheet and annotate which seqs originate from there 
# For me this easy easy / quick, but if this is not ideal for your data you can generate a script to pull the origin as well 
###

########################################################################
# STEP 3: 
# Running the PairwiseSoftAlignment tool 
# installation instructions: 
########################################################################

python /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/sa_interface.py /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0013/cat_GVOGm0013_cleaned.faa /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0013/GVOGm0013_output.tsv

python /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/sa_interface.py /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0022/cat_GVOGm0022_cleaned.faa /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0022/GVOGm0022_output.tsv

python /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/sa_interface.py /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0023/cat_GVOGm0023_cleaned.faa /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0023/GVOGm0023_output.tsv

python /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/sa_interface.py /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0461/cat_GVOGm0461_cleaned.faa /Users/syrenawhitner/Desktop/Will_tool/PairwiseSoftAlignments-main/Round2/GVOGm0461/GVOGm0461_output.tsv

### 
# NOTE : the output TSV will go into the R script "PCoA.R"


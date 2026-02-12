########################################################################
#
# Predicting and annotating viral / fungal proteins
# Viral Marker Classification  
# This notebook utilizes several GV reference Hidden Markov Models (hmms) sets, as referred to: 
# "Core viral" - https://github.com/faylward/ncldv_markersearch/blob/master/hmm/NCLDV.hmm (Moniruzzaman, 
# Martinez-Gutierrez, et al., 2020)
# "Common viral" - https://github.com/faylward/ncldv_markersearch/blob/master/hmm/gvogs.common.hmm 
# (Moniruzzaman, Martinez-Gutierrez, et al., 2020)
# "Full NCLDV" - https://faylward.github.io/GVDB/ (Aylward et al., 2021)
# Syrena Whitner
# December, 2025
#
########################################################################
# Step 1: 
# Viral protein/ORF prediction using Prodigal & Prodigal-GV for all candidate contigs 
# Candidate contigs generated in Genomad_CandidateContigs notebook
########################################################################
# with Prodigal (v2.6.3)
########################################################################

# running prodigal 
# Directories
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/sequences"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/prodigal_output"

for file in "$indir"/*.fasta; do
    base=$(basename "$file" .fasta)

    echo "Running Prodigal on $base..."

    prodigal \
        -i "$file" \
        -o "${outdir}/${base}.gbk" \
        -a "${outdir}/${base}.faa" \
        -d "${outdir}/${base}.fna" \
        -f gbk \
        -p meta
done


########################################################################
# With Prodigal-GV (v2.11.0)
########################################################################
# running prodigal-gv
# Directories
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/sequences"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/prodigalGV_output"

for file in "$indir"/*.fasta; do
    base=$(basename "$file" .fasta)

    echo "Running Prodigal on $base..."

    prodigal-gv \
        -i "$file" \
        -o "${outdir}/${base}.gbk" \
        -a "${outdir}/${base}.faa" \
        -d "${outdir}/${base}.fna" \
        -f gbk \
        -p meta
done


########################################################################
# Step 2: 
# Predicted eukaryotic gene models (CDS, introns, splicing) using AUGUSTUS (v3.5.0) pretrained with Aspergillus nidulans
########################################################################

#!/bin/bash
#SBATCH --job-name=augustus_predictions
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=03-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --array=1-10   
#SBATCH --cpus-per-task=12
#SBATCH --mem=400G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

ml lang/Anaconda3/2024.02-1
module load tools/parallel/20220722-GCCcore-11.3.0
source activate funannotate

# Directories
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/phylogeny/funannotate_scaffolds/masked_reads"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/augustus_output"
mkdir -p "$outdir" logs

# make a list of input files
genomes=($(ls "$indir"/*.fasta))
genome=${genomes[$SLURM_ARRAY_TASK_ID-1]}
base=$(basename "$genome" .fasta)

echo "[$(date)] Starting AUGUSTUS on $base"

augustus --species=aspergillus_nidulans \
         --softmasking=1 \
         --gff3=on \
         "$genome" > "$outdir/${base}.augustus.gff3" 2> "$outdir/${base}.log"

echo "[$(date)] Finished $base"

# also ran on just the subset of contigs that had "viral reads"
# Directories
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/sequences"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/augustus_output"
mkdir -p "$outdir" logs

for genome in "$indir"/*.fasta; do
    base=$(basename "$genome" .fasta)
    echo "[$(date)] Starting AUGUSTUS on $base"

    augustus --species=aspergillus_nidulans \
             --softmasking=1 \
             --gff3=on \
             "$genome" > "$outdir/${base}.augustus.gff3" \
             2> "$outdir/${base}.log"

    echo "[$(date)] Finished $base"
done


########################################################################
# Step 3: 
#  Viral marker identification
# Searched predicted viral ORFs against NCLDV-relevant HMM databases using HMMER3 (v3.4)
########################################################################
# Queried the "Core" NCLDV hallmark marker set 
# https://github.com/faylward/ncldv_markersearch/blob/master/hmm/NCLDV.hmm 
# (Moniruzzaman, Martinez-Gutierrez, et al., 2020)
########################################################################

for hmm in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/markersearch_set/NCLDV_hmm_set/*.hmm; do
  hmmname=$(basename "$hmm" .hmm)

  for faa in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/syrena_seqs/prodigal_output/*.faa; do
    faname=$(basename "$faa" .faa)

    outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/syrena_seqs/hmm_output_ncldv"
    mkdir -p "$outdir"

    tblout="${outdir}/${faname}_${hmmname}.tbl"
    domout="${outdir}/${faname}_${hmmname}.domtbl"

    echo "Running hmmsearch on $faname with $hmmname..."
    hmmsearch --cpu 12 \
              --tblout "$tblout" \
              --domtblout "$domout" \
              "$hmm" "$faa" > "${outdir}/${faname}_${hmmname}.log"
  done
done

########################################################################
# Queried the "Common" occurring GVOG marker set
# https://github.com/faylward/ncldv_markersearch/blob/master/hmm/gvogs.common.hmm 
# (Moniruzzaman, Martinez-Gutierrez, et al., 2020)
########################################################################

for hmm in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/markersearch_set/gvogs_common_hmm/*.hmm; do
  hmmname=$(basename "$hmm" .hmm)

  for faa in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/prodigal_output/*.faa; do
    faname=$(basename "$faa" .faa)

    outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/hmm_output_common"
    mkdir -p "$outdir"

    tblout="${outdir}/${faname}_${hmmname}.tbl"
    domout="${outdir}/${faname}_${hmmname}.domtbl"

    echo "Running hmmsearch on $faname with $hmmname..."
    hmmsearch --cpu 8 \
              --tblout "$tblout" \
              --domtblout "$domout" \
              "$hmm" "$faa" > "${outdir}/${faname}_${hmmname}.log"
  done
done

########################################################################
# Queried the "full" GVOG database
# https://faylward.github.io/GVDB/ 
# (Aylward et al., 2021)
########################################################################

for hmm in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/GVOGs/*.hmm; do
  hmmname=$(basename "$hmm" .hmm)

  for faa in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/prodigal_output/*.faa; do
    faname=$(basename "$faa" .faa)

    out="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/hmm_output_fullset/${faname}_${hmmname}.tbl"
    echo "Running hmmsearch on $faname with $hmmname..."
    hmmsearch --tblout "$out" "$hmm" "$faa"
  done
done


########################################################################
# Parsed and filtered dom tables 
# Extracted matching protein sequences using seqkit (v2.10.1) -> 
# this was only done from the "common" hmms because these seemed the most sensitive 
# these were what was used for phylogeny building downstream 
########################################################################

# making a top hits table and extracting these 
# Directories
RES_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/hmm_output_common"
OUT_SUMMARY="$RES_DIR/best_hits_summary.tsv"

echo "Parsing HMM results from: $RES_DIR"

# Initialize summary header
echo -e "Genome\tHMM\tHit\tEvalue\tScore" > "$OUT_SUMMARY"

# Loop through all .tbl files
for tbl in "$RES_DIR"/*.tbl; do
  genome=$(basename "$tbl" | sed 's/_GVOGm[0-9]*\.trim\.tbl//')
  hmm=$(basename "$tbl" | sed -E 's/.*_(GVOGm[0-9]+)\.trim\.tbl/\1/')

  # grab the first non-comment line (best hit)
  line=$(grep -v "^#" "$tbl" | head -n 1 || true)
  if [ -n "$line" ]; then
    hit=$(echo "$line" | awk '{print $1}')
    evalue=$(echo "$line" | awk '{print $5}')
    score=$(echo "$line" | awk '{print $6}')
  else
    hit="NA"; evalue="NA"; score="NA"
  fi

  echo -e "$genome\t$hmm\t$hit\t$evalue\t$score" >> "$OUT_SUMMARY"
done

echo "Wrote raw summary: $OUT_SUMMARY"

##### extracting the top hit from each original protein file 
# Directories 
BASE="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL"
FAA_DIR="$BASE/prodigal_output"
TABLE="$BASE/hmm_output_common/best_hits_summary.tsv"
OUT_DIR="$BASE/hmm_extracts"

mkdir -p "$OUT_DIR"

echo "Using input table: $TABLE"
echo "Protein directory: $FAA_DIR"
echo "Output directory: $OUT_DIR"
echo

tail -n +2 "$TABLE" | while IFS=$'\t' read -r genome hmm hit evalue score; do
  # Skip NA or empty hits
  if [[ "$hit" == "NA" || -z "$hit" ]]; then
    continue
  fi

  # Search for the file (case-insensitive, ends in .faa)
  faa_file=$(find "$FAA_DIR" -type f -iname "${genome}*.faa" | head -n 1)

  out_file="$OUT_DIR/${genome}_${hmm}_besthit.faa"

  if [[ -n "$faa_file" && -f "$faa_file" ]]; then
    echo "Extracting $hit from $(basename "$faa_file") → $(basename "$out_file")"
    seqkit grep -p "$hit" "$faa_file" > "$out_file" 2>/dev/null || \
      echo "Warning: Could not find $hit in $(basename "$faa_file")"
  else
    echo "No matching FASTA found for genome: $genome"
  fi
done

########################################################################
# Step 4: MANUAL CURATION - This completely depends on your end goal so is not "plug and chug"
# 
# Check outputs manually to see which isolates pass thresholds for most "core" NCLDV markers, 
# Constructed a GVOG presence/absence matrix in R
# Filtered HMM hits using an e-value threshold of 1×10⁻³
# Identified GVOGs most frequently detected across isolates
# top 8 most occuring gvogs, and then identified genomes in which any combination of these occur at least 4 times 
# GVOGm 0013, 0020, 0022, 0023, 0054, 0214, 0461, and 1574. 
# from these genomes 
"1002337_sa182576_combined"
"1002442_sa182533_combined"
"1002610_sa182534_combined"
"1002970_sa182557_combined"
"1003183_sa150016_combined"
"1003189_sa149994_combined"
"1003356_sa149980_combined"
"UH_1000136_sa123987_combined"
"UH_1000339_sa123914_combined"
"UH_1000590_sa144727_combined"
"UH_1001061_sa129894_combined"
"UH_1001187_sa123975_combined"
"UH_1001220_sa145563_combined"
"UH_1001271_sa144725_combined"
"UH_1001376_sa129866_combined"
"UH_1001502_sa132352_combined"
"UH_1001523_sa129663_combined"
"UH_1001535_sa132328_combined"
"UH_1001556_sa132306_combined"
"UH_1001604_sa129678_combined"
"UH_1002031_sa145543_combined"
########################################################################
# now putting them into their own directory 
########################################################################

# Directories
src="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/hmm_extracts"
dest="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/top8_gvog_extracts"

# Create destination directory if it doesn't exist
mkdir -p "$dest"

# Define target GVOGs
gvogs=("GVOGm0013" "GVOGm0020" "GVOGm0022" "GVOGm0023" "GVOGm0054" "GVOGm0214" "GVOGm0461" "GVOGm1574")

# Define your target genomes
genomes=(
"1002337_sa182576_combined"
"1002442_sa182533_combined"
"1002610_sa182534_combined"
"1002970_sa182557_combined"
"1003183_sa150016_combined"
"1003189_sa149994_combined"
"1003356_sa149980_combined"
"UH_1000136_sa123987_combined"
"UH_1000339_sa123914_combined"
"UH_1000590_sa144727_combined"
"UH_1001061_sa129894_combined"
"UH_1001187_sa123975_combined"
"UH_1001220_sa145563_combined"
"UH_1001271_sa144725_combined"
"UH_1001376_sa129866_combined"
"UH_1001502_sa132352_combined"
"UH_1001523_sa129663_combined"
"UH_1001535_sa132328_combined"
"UH_1001556_sa132306_combined"
"UH_1001604_sa129678_combined"
"UH_1002031_sa145543_combined"
)

# Copy matching files
for genome in "${genomes[@]}"; do
  for gvog in "${gvogs[@]}"; do
    file="${src}/${genome}_${gvog}_besthit.faa"
    if [[ -f "$file" ]]; then
      echo "📂 Copying: $(basename "$file")"
      cp "$file" "$dest/"
    else
      echo "Missing: $file"
    fi
  done
done

echo "Done! All available best hits copied to: $dest"

# renaming headers 
for f in /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/top8_gvog_extracts/*.faa; do
  base=$(basename "$f" _combined_*.faa)       # strip after "_combined_"
  genome="${base%%_combined*}"                # remove everything after "_combined"
  echo "Renaming headers in $f → >$genome"
  sed -i "s/^>.*/>${genome}/" "$f"
done

# combine markers 
i.e - 
cat *_GVOGm1574_* > ./cat_markers/GVOGm1574_cat.faa


########################################################################
# Step 5: 
# Viral marker identification from other reference NCLDVs
# Here I am pulling from other known NCLDVs that infect fungi
# Myers et al. 2025 
# "Mycodnoviridae" https://deepblue.lib.umich.edu/data/concern/data_sets/9880vr84k
#
# also pulling from ICTV reference NCLDVs 
# https://faylward.github.io/GVDB/
#
# prior to this analysis, both datasets were queried using the "common" hmm marker set with prodigal as in Step #1
########################################################################
# Note I am only searching for specific hmms after curating set from previous step 
#### for umich files 
HMM_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/markersearch_set/gvogs_common_hmm"
FAA_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/umich_gv/prodigal_output"
OUT_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs/umich_hmm_output"

mkdir -p "$OUT_DIR"

# Target HMMs 
hmms=("GVOGm0020" "GVOGm0022" "GVOGm0054" "GVOGm0214" "GVOGm1574")

# Loop
for hmm in "${hmms[@]}"; do
  hmmfile="$HMM_DIR/${hmm}.trim.hmm"
  if [[ ! -f "$hmmfile" ]]; then
    echo "Missing HMM file: $hmmfile"
    continue
  fi

  for faa in "$FAA_DIR"/*.faa; do
    faname=$(basename "$faa" .faa)
    outbase="${OUT_DIR}/${faname}_${hmm}"

    echo "Running hmmsearch on $faname with $hmm ..."
    hmmsearch --cpu 8 \
              --tblout "${outbase}.tbl" \
              --domtblout "${outbase}.domtbl" \
              "$hmmfile" "$faa" > "${outbase}.log"
  done
done

echo "Done! Results saved to: $OUT_DIR"

# Directories
RES_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs/umich_hmm_output"
OUT_SUMMARY="$RES_DIR/best_hits_summary.tsv"

echo "Parsing HMM results from: $RES_DIR"

# Initialize summary header 
echo -e "Genome\tHMM\tHit\tEvalue\tScore" > "$OUT_SUMMARY"

# Loop through all .tbl files 
for tbl in "$RES_DIR"/*.tbl; do
  genome=$(basename "$tbl" | sed 's/_GVOGm[0-9]*\.trim\.tbl//')
  hmm=$(basename "$tbl" | sed -E 's/.*_(GVOGm[0-9]+)\.trim\.tbl/\1/')

  # grab the first non-comment line (best hit)
  line=$(grep -v "^#" "$tbl" | head -n 1 || true)
  if [ -n "$line" ]; then
    hit=$(echo "$line" | awk '{print $1}')
    evalue=$(echo "$line" | awk '{print $5}')
    score=$(echo "$line" | awk '{print $6}')
  else
    hit="NA"; evalue="NA"; score="NA"
  fi

  echo -e "$genome\t$hmm\t$hit\t$evalue\t$score" >> "$OUT_SUMMARY"
done

echo "Wrote raw summary: $OUT_SUMMARY"


# for some reason the naming is weird so need to fixx the .tsv file 
awk 'BEGIN{OFS="\t"}
NR==1 { print; next }  # keep header
{
  # sanitize CRs (just in case)
  gsub(/\r/,"")

  # 1) Clean Genome: remove trailing _GVOGm####(.trim)?.tbl  OR lone .tbl
  genome = $1
  sub(/_GVOGm[0-9]+(\.trim)?\.tbl$/, "", genome)
  sub(/\.tbl$/, "", genome)            # fallback if needed
  $1 = genome

  # 2) Clean HMM: keep only GVOGm####
  if (match($2, /(GVOGm[0-9]+)/, m)) $2 = m[1]

  print
}' best_hits_summary.tsv > filtered_best_hits_summary.tsv

##### extracting the top hit from each original protein file 
# Directories
BASE="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs"
FAA_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/umich_gv/prodigal_output"
TABLE="$BASE/umich_hmm_output/filtered_best_hits_summary.tsv"
OUT_DIR="$BASE/umich_hmm_extracts"

mkdir -p "$OUT_DIR"

tail -n +2 "$TABLE" | while IFS=$'\t' read -r genome hmm hit evalue score; do
  # Skip empty or NA hits
  [[ "$hit" == "NA" || -z "$hit" ]] && continue

  faa_file=$(find "$FAA_DIR" -type f -iname "${genome}.faa" | head -n 1)
  out_file="$OUT_DIR/${genome}_${hmm}_besthit.faa"

  if [[ -n "$faa_file" && -f "$faa_file" ]]; then
    echo "▶️ Extracting $hit from $(basename "$faa_file") → $(basename "$out_file")"

    # Extract record whose header begins with the hit (before first space)
    awk -v id="$hit" '
      /^>/ {
        if (seq) {
          print seq
        }
        seq = ""
        header = $0
        split($0, a, " ")
        if (a[1] == ">"id) {
          print header
          getline
          while ($0 !~ /^>/ && NF) {
            print
            if (!getline) break
          }
          exit
        }
      }
    ' "$faa_file" > "$out_file"

    nheads=$(grep -c "^>" "$out_file" || true)
    if [[ "$nheads" -ne 1 ]]; then
      echo "Warning: expected 1 record for $hit, got $nheads in $(basename "$out_file")"
    fi
  else
    echo "No matching FASTA found for genome: $genome"
  fi
done

echo "Extraction complete!"
echo "Output directory: $OUT_DIR"


# rename headers to match parent file 
DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs/umich_hmm_extracts"

echo "Renaming FASTA headers (keeping only genome name) in: $DIR"
echo

for f in "$DIR"/*_besthit.faa; do
  base=$(basename "$f" _besthit.faa)
  # remove the last underscore + token (the HMM name)
  genome=$(echo "$base" | sed -E 's/_[^_]+$//')
  tmp="${f}.tmp"

  echo "Processing: $(basename "$f") → header = ${genome}"
  awk -v newhdr=">${genome}" '
    /^>/ {print newhdr; next}
    {print}
  ' "$f" > "$tmp" && mv "$tmp" "$f"
done

echo
echo "Header renaming complete!"
echo "All headers now use genome names only."


cp *_cat.faa /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs/cat_markers

cat *_GVOGm1574_* > /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/umich_seqs/cat_markers/GVOGm1574_cat.faa

########################################################################
# now for ictv references 
########################################################################

# Directories 
HMM_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/markersearch_set/gvogs_common_hmm"
FAA_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/GVDB_genomes/GVDB_prot"
OUT_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_hmm_output"

mkdir -p "$OUT_DIR"

# Target HMMs 
hmms=("GVOGm0020" "GVOGm0214" "GVOGm1574")

# Loop 
for hmm in "${hmms[@]}"; do
  hmmfile="$HMM_DIR/${hmm}.trim.hmm"
  if [[ ! -f "$hmmfile" ]]; then
    echo "Missing HMM file: $hmmfile"
    continue
  fi

  for faa in "$FAA_DIR"/*.faa; do
    faname=$(basename "$faa" .faa)
    outbase="${OUT_DIR}/${faname}_${hmm}"

    echo "Running hmmsearch on $faname with $hmm ..."
    hmmsearch --cpu 8 \
              --tblout "${outbase}.tbl" \
              --domtblout "${outbase}.domtbl" \
              "$hmmfile" "$faa" > "${outbase}.log"
  done
done

echo "Done! Results saved to: $OUT_DIR"

# Directories 
RES_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_hmm_output"
OUT_SUMMARY="$RES_DIR/best_hits_summary.tsv"

echo "Parsing HMM results from: $RES_DIR"

# Initialize summary header
echo -e "Genome\tHMM\tHit\tEvalue\tScore" > "$OUT_SUMMARY"

# Loop through all .tbl files
for tbl in "$RES_DIR"/*.tbl; do
  genome=$(basename "$tbl" | sed 's/_GVOGm[0-9]*\.trim\.tbl//')
  hmm=$(basename "$tbl" | sed -E 's/.*_(GVOGm[0-9]+)\.trim\.tbl/\1/')

  # grab the first non-comment line (best hit)
  line=$(grep -v "^#" "$tbl" | head -n 1 || true)
  if [ -n "$line" ]; then
    hit=$(echo "$line" | awk '{print $1}')
    evalue=$(echo "$line" | awk '{print $5}')
    score=$(echo "$line" | awk '{print $6}')
  else
    hit="NA"; evalue="NA"; score="NA"
  fi

  echo -e "$genome\t$hmm\t$hit\t$evalue\t$score" >> "$OUT_SUMMARY"
done

echo "Wrote raw summary: $OUT_SUMMARY"

# for some reason the naming is messed up so need to fixx the .tsv file 
awk 'BEGIN{OFS="\t"}
NR==1 { print; next }  # keep header
{
  # sanitize CRs (just in case)
  gsub(/\r/,"")

  # 1) Clean Genome: remove trailing _GVOGm####(.trim)?.tbl  OR lone .tbl
  genome = $1
  sub(/_GVOGm[0-9]+(\.trim)?\.tbl$/, "", genome)
  sub(/\.tbl$/, "", genome)            # fallback if needed
  $1 = genome

  # 2) Clean HMM: keep only GVOGm####
  if (match($2, /(GVOGm[0-9]+)/, m)) $2 = m[1]

  print
}' best_hits_summary.tsv > filtered_best_hits_summary.tsv

##### extracting the top hit from each original protein file 
# Directories
BASE="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_seqs"
FAA_DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/NCLDV_files/GVDB_genomes/GVDB_prot"
TABLE="$BASE/ictv_hmm_output/filtered_best_hits_summary.tsv"
OUT_DIR="$BASE/ictv_hmm_extracts"

mkdir -p "$OUT_DIR"


tail -n +2 "$TABLE" | while IFS=$'\t' read -r genome hmm hit evalue score; do
  [[ "$hit" == "NA" || -z "$hit" ]] && continue

  hit=$(echo "$hit" | tr -d '\r')

  faa_file=$(find "$FAA_DIR" -type f \( -iname "${genome}*.faa" -o -iname "${genome}*.fa" \) | head -n 1)
  out_file="$OUT_DIR/${genome}_${hmm}_besthit.faa"

  if [[ -n "$faa_file" && -f "$faa_file" ]]; then
    echo "Extracting $hit from $(basename "$faa_file") → $(basename "$out_file")"

    # Escape any regex-sensitive characters in the hit ID
    safe_hit=$(printf '%s\n' "$hit" | sed 's/[][^$.*/]/\\&/g')

    # Match start (^) + literal hit ID, stop at space or #
    seqkit grep -r -p "^${safe_hit}[[:space:]#]" "$faa_file" \
      | awk 'BEGIN{RS=">"; ORS=""} NR==2{print ">"$0}' > "$out_file"

    nheads=$(grep -c "^>" "$out_file" || true)
    if [[ "$nheads" -ne 1 ]]; then
      echo "Warning: expected 1 record for $hit, got $nheads in $(basename "$out_file")"
    fi
  else
    echo "No matching protein FASTA found for genome: $genome"
  fi
done

echo "Extraction complete!"
echo "Output directory: $OUT_DIR"


# rename headers to match parent file 
DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_seqs/cat_markers"

echo "🔧 Cleaning FASTA headers (removing everything after first '|')..."
echo

for f in "$DIR"/*.faa; do
  echo "Fixing $(basename "$f")"
  tmp="${f}.tmp"

  # Correct: use extended regex parentheses without backslashes
  sed -E 's/^(>[^|]*)\|.*$/\1/' "$f" > "$tmp" && mv "$tmp" "$f"
done

echo
echo "Headers cleaned!"
echo "All text after '|' removed from FASTA headers."

DIR="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_seqs/cat_markers"

for f in "$DIR"/GVOGm0022.copy1.faa "$DIR"/GVOGm0054.copy1.faa; do
  if [[ -f "$f" ]]; then
    echo "Fixing $(basename "$f")"
    awk 'BEGIN{OFS=""} 
         /^>/ { print ">"$NF; next }  # header: set to ">" + last field
              { print }               # sequence lines: print as-is
        ' "$f" > "$f.tmp" && mv "$f.tmp" "$f"
  else
    echo "Missing: $f"
  fi
done

# now need to cat and move over the other markers 
cat *_GVOGm1574_* > /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/ictv_seqs/cat_markers/GVOGm1574_cat.faa

########################################################################
# Step 6:
# Just putting everything in one place 
# at this point, I should have a master file for each candidate GVOG that has 
# my sequences, umich references, and from ICTV as well 
# these will be used in phylogenetic construction 
########################################################################

# now I have the directory all_seqs which is where I will put everything! 
# Directories 
BASE="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL"

DIR1="$BASE/syrena_seqs/top8_gvog_extracts/cat_markers" #took seqs from my own previous steps, put everything into syrena_seqs just to keep seaprate 
DIR2="$BASE/ictv_seqs/cat_markers"
DIR3="$BASE/umich_seqs/cat_markers"

OUT_DIR="$BASE/all_seqs/cat_markers"
mkdir -p "$OUT_DIR"

echo "Combining marker sequences from:"
echo "   - $DIR1"
echo "   - $DIR2"
echo "   - $DIR3"
echo "Output: $OUT_DIR"
echo

# Loop over all marker files in the first directory
for marker in "$DIR1"/*.faa; do
  file=$(basename "$marker")                # e.g. GVOGm0013_cat.faa
  out_file="$OUT_DIR/$file"

  echo "Combining sequences for: $file"

  # Concatenate the same file from each directory if it exists
  for d in "$DIR1" "$DIR2" "$DIR3"; do
    if [[ -f "$d/$file" ]]; then
      cat "$d/$file" >> "$out_file"
    fi
  done

  echo "Done: $out_file"
  echo
done

echo "Combination complete!"
echo "One combined file per GVOG marker is in: $OUT_DIR"
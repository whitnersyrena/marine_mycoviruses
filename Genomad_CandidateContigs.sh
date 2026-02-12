########################################################################
#
# Initial screening for NCLDV-like sequences
# This notebook covers the initial screening of fungal genomes for viral 
# sequences using GENOMAD (v1.11.2), and then additional contig curation based on 
# values such as gv_marker, taxonomic classification, etc 
# Syrena Whitner
# December, 2025
#
########################################################################
# STEP 1: 
# Screened assembled contigs using geNomad (end-to-end workflow)
########################################################################

    #!/bin/bash
    #SBATCH --job-name=20230919_trial_sortmerna 
    #SBATCH --partition=shared

    ## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
    #SBATCH --time=00-08:00:00 ## time format is DD-HH:MM:SS
    #SBATCH --nodes=1
    #SBATCH --cpus-per-task=4
    #SBATCH --mem=32G ## max amount of memory per node you require

    #SBATCH --error=job%A.err ## %A - filled with jobid
    #SBATCH --output=job%A.out ## %A - filled with jobid
    #SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
    #SBATCH --mail-user=swhitner@hawaii.edu

    ml lang/Anaconda3/2022.05
    source activate genomad 
    cd /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset

    for sample in $(cat sample_names.txt); do genomad end-to-end ${sample}.fasta output trial_output; done 

########################################################################
# Step 2: 
# Identified contigs classified as Nucleocytoviricota or Unclassified
# Retained contigs with GV marker frequency >20%
# Concatenated retained candidate contigs per fungal isolate
########################################################################

# copying all the features to one folder 
find /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset \
  -type f -path "*/*_marker_classification/*_features.tsv" \
  -exec cp {} /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/features/ \;

# now pulling the summaries 
find /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset \
  -type f -path "*/*_summary/*_summary.tsv" \
  -exec cp {} /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/summaries/ \;

# making the .txt files for candidate sequences 
# folder with  *_features.tsv files
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/features"

for f in "$indir"/*_features.tsv; do
  out="${f%_features.tsv}_sequences.txt"
  awk -F'\t' '
    FNR==1 {
      # find column indices by exact header match
      gv_idx = seq_idx = 0
      for (i=1; i<=NF; i++) {
        if ($i == "gv_marker_freq") gv_idx = i
        if ($i == "seq_name")      seq_idx = i
      }
      next
    }
    gv_idx && seq_idx && ($gv_idx + 0 >= 0.20) {
      print $seq_idx
    }
  ' "$f" > "$out"
done

# doing the same for the summaries 
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/summaries"


for f in "$indir"/*_summary.tsv; do
  out="${f%_summary.tsv}_sequences.txt"
  awk -F'\t' '
    FNR==1 {
      # find column indices
      tax_idx = seq_idx = 0
      for (i=1; i<=NF; i++) {
        if ($i == "taxonomy") tax_idx = i
        if ($i == "seq_name") seq_idx = i
      }
      next
    }
    tax_idx && seq_idx && ($tax_idx == "Unclassified" || $tax_idx ~ /^Viruses;Varidnaviria;Bamfordvirae;Nucleocytoviricota/) {
      print $seq_idx
    }
  ' "$f" > "$out"
done

# combining them into one single sequence list 
features_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/features"
summaries_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/summaries"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/sequence_lists"

mkdir -p "$outdir"

# Normalize summary filenames by removing "_virus" before merging
for f in "$summaries_dir"/*_virus_sequences.txt; do
  base=$(basename "$f" _virus_sequences.txt)
  newname="${summaries_dir}/${base}_sequences.txt"
  mv "$f" "$newname"
done

# Merge matching files from both directories and remove duplicates
for f in "$features_dir"/*_sequences.txt; do
  base=$(basename "$f")
  summary_file="$summaries_dir/$base"
  outfile="$outdir/$base"

  # If a corresponding summary file exists, combine and deduplicate
  if [[ -f "$summary_file" ]]; then
    cat "$f" "$summary_file" | sort -u > "$outfile"
  else
    # If not, just copy the features file
    sort -u "$f" > "$outfile"
  fi
done

# Handle any summary files that didn’t have a matching feature file
for f in "$summaries_dir"/*_sequences.txt; do
  base=$(basename "$f")
  outfile="$outdir/$base"
  [[ -f "$outfile" ]] || sort -u "$f" > "$outfile"
done

# now puttling all these contigs from each parent genome, putting into directory ./genomad_output_subset/candidate_conitgs
seq_list_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/sequence_lists"
genome_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/hexagon_genomes/genomes/for_viruses"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/candidate_contigs"

mkdir -p "$outdir"

# Loop through every *_sequences.txt list
for list in "$seq_list_dir"/*_sequences.txt; do
  base=$(basename "$list" _sequences.txt)
  genome="$genome_dir/${base}.fasta"
  outfile="$outdir/${base}_subset.fasta"

  if [[ ! -f "$genome" ]]; then
    echo "Genome not found for $base — skipping"
    continue
  fi

  echo "Extracting contigs for $base ..."

  # Build an awk pattern file (quoted list of contig IDs)
  awk '{print "^>"$1"$"}' "$list" > "${list%.txt}_patterns.tmp"

  # Extract matching contigs (+ sequences) from the genome fasta
  awk -v patfile="${list%.txt}_patterns.tmp" '
    BEGIN {
      while ((getline p < patfile) > 0) pats[p] = 1
    }
    /^>/ {
      keep = 0
      for (p in pats) if ($0 ~ p) { keep = 1; break }
    }
    keep
  ' "$genome" > "$outfile"

  rm -f "${list%.txt}_patterns.tmp"
done

echo "Done! Extracted contigs are in: $outdir"

# combined my candidate sequences with the highest confidence predictions from genomad as well 

# Directories
# dir1 is the output frome genomad 
dir1="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/virus_genomes"
dir2="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/candidate_contigs/sequences"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/sequences"

mkdir -p "$outdir"

# Loop through all *_virus.fna files
for virus in "$dir1"/*_virus.fna; do
    base=$(basename "$virus" _virus.fna)
    subset="$dir2/${base}_subset.fasta"
    outfile="$outdir/${base}_combined.fasta"

    if [[ -f "$subset" ]]; then
        echo "🧬 Combining $base"
        cat "$subset" "$virus" > "$outfile"
    else
        echo "No matching subset file for $base"
    fi
done

echo "Done! Combined FASTAs written to $outdir"
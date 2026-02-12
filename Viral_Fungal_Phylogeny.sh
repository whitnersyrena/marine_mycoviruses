########################################################################
#
# Viral phylogenetic construction
# Fungal phylogenetic construction
# 
# Syrena Whitner
# December, 2025
#
########################################################################

########################################################################
# Step #1:
# Constructing VIRAL phylogeny using the following four genes:
# GVOGm00 13, 22, 23, and 461 - will need at least 2/4 of the genes present 
# removing the seqs that only have one of the markers present 
########################################################################

# directories 
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/all_seqs/cat_markers"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/cat_markers"

mkdir -p "$outdir"

# list of .faa files
files=(
  GVOGm0013_cat.faa
  GVOGm0022_cat.faa
  GVOGm0023_cat.faa
  GVOGm0461_cat.faa
)

# equence headers to remove (without ">") 
remove_ids=(
  "1002337_sa182576"
  "1002970_sa182557"
  "UH_1001604_sa129678"
)

# Loop through each file 
for f in "${files[@]}"; do
  infile="${indir}/${f}"
  outfile="${outdir}/${f}"

  echo "Cleaning $f ..."
  awk -v RS='>' -v ORS='' -v ids="$(IFS='|'; echo "${remove_ids[*]}")" '
    BEGIN {
      split(ids, arr, "|")
      for (i in arr) bad[arr[i]] = 1
    }
    NR > 1 {
      header = $1
      split(header, parts, /[ \t\n]/)
      id = parts[1]
      if (!(id in bad)) print ">" $0
    }' "$infile" > "$outfile"

  echo "Wrote cleaned file: $outfile"
done

echo "Done. Sequences removed: ${remove_ids[*]}"

########################################################################
# Also removing some sequences from reference NCLDV genomes from the ICTV giant virus database
# just to speed up phylogeny construction
########################################################################

# directories 
indir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/cat_markers"
outdir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/cat_markers_imt.removed"
remove_list="$indir/remove_list.txt"

mkdir -p "$outdir"

# Check that remove_list exists 
if [[ ! -f "$remove_list" ]]; then
  echo "Error: remove_list.txt not found in $indir"
  exit 1
fi

# Filtering loop 
for file in "$indir"/*.faa; do
  base=$(basename "$file")
  echo "🔹 Filtering $base ..."

  # AWK logic:
  # 1. Read the list of headers (starting with >)
  # 2. Remove any sequence where the header exactly matches one in the list
  awk -v rmfile="$remove_list" '
    BEGIN {
      while ((getline line < rmfile) > 0) {
        gsub(/\r/, "", line)        # remove Windows line endings if present
        if (line ~ /^>/) bad[line] = 1
      }
    }
    /^>/ {
      keep = !($0 in bad)
    }
    keep
  ' "$file" > "$outdir/$base"

  echo "Saved filtered file: $outdir/$base"
done

echo "All files processed. Filtered versions are in:"
echo "   $outdir"

########################################################################
# Step 2: 
# Aligned individual GVOGs using MAFFT (v7.5)
# Trimmed alignments using TrimAl (v1.5) with automated settings
# Concatenated trimmed alignments into a supermatrix using AMAS (v1.2)
# Generated partition files corresponding to individual GVOGs
# Inferred a concatenated viral phylogeny using IQ-TREE2 (v2.4.0)
# Performed automated model selection and partition merging
# Assessed branch support with ultrafast bootstrap and SH-aLRT tests
########################################################################

#!/bin/bash
#SBATCH --job-name=allvir_4gvog_full.sh
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=07-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=600G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

# Activate environment and load modules 
ml lang/Anaconda3/2024.02-1
source activate viral_phylogeny
module load bio/MAFFT/7.505-GCC-11.3.0-with-extensions
cd /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5

# Directories 
base="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/cat_markers_imt.removed"
aln_dir="$base/aln"
trim_dir="$base/aln_trimmed"
matrix_dir="$base/matrix"

mkdir -p "$aln_dir" "$trim_dir" "$matrix_dir"

# 1. Align with MAFFT -> GVOGmXXXX_cat.aln.faa 
echo "Aligning sequences with MAFFT..."
for gvog in "$base"/GVOGm*.faa; do
  stem=$(basename "$gvog" .faa)                    # e.g. GVOGm0022_cat
  mafft --auto "$gvog" > "$aln_dir/${stem}.aln.faa"
done
echo "MAFFT step complete."

# 2. Trim with trimAl -> GVOGmXXXX_cat.aln.trim.faa 
echo "Trimming alignments..."
for aln in "$aln_dir"/*.aln.faa; do
  # Take '...aln.faa', drop the final '.faa', then append '.trim.faa'
  # Result: '...aln.trim.faa'
  out_name="$(basename "${aln%.faa}.trim.faa")"
  trimal -in "$aln" -out "$trim_dir/$out_name" -automated1
done
echo "trimAl step complete."

# 3. Concatenate with AMAS 
echo "Concatenating alignments..."
AMAS concat \
  -f fasta \
  -d aa \
  -i "$trim_dir"/*.aln.trim.faa \
  -p "$matrix_dir/gvog_partitions.txt" \
  -t "$matrix_dir/gvog_supermatrix.faa" \
  -y raxml
echo "AMAS concatenation complete."

# 4. build tree
echo "Building ML tree with IQ-TREE3..."
iqtree3 \
  -s "$matrix_dir/gvog_supermatrix.faa" \
  -p "$matrix_dir/gvog_partitions.txt" \
  -T 16 \
  -m MFP+MERGE \
  -B 1000 \
  -alrt 1000 \
  --prefix "$matrix_dir/myseqs_robust_tree" \
  -v

echo "Pipeline finished successfully! - tree is ready :) "


########################################################################
# Step 3: 
# Fungal phylogeny construction 
# Selected fungal isolates meeting viral hallmark inclusion criteria
# Identified conserved fungal marker genes using BUSCO Fungi_odb10
# Generated gene-wise protein alignments using Phyling (v2.3.0)
# Constructed a concatenated fungal phylogeny using IQ-TREE2 via Phyling
# Performed partitioned model selection and branch support estimation
########################################################################

# starting with the 4gvog subset

# Define genome list --> these are the same genomes I use for my viral analysis
genomes=(
  "1002442_sa182533"  "1002610_sa182534" "1003183_sa150016" "1003189_sa149994"
  "1003356_sa149980"  "UH_1000136_sa123987" "UH_1000339_sa123914" "UH_1000590_sa144727"
  "UH_1001061_sa129894" "UH_1001187_sa123975" "UH_1001220_sa145563" "UH_1001271_sa144725"
  "UH_1001376_sa129866" "UH_1001502_sa132352" "UH_1001523_sa129663" "UH_1001535_sa132328"
  "UH_1001556_sa132306" "UH_1002031_sa145543"
)

# directories 
src="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/phylogeny/prodigal_results/fungal_proteins"
dest="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal/fun_proteins"

# Make sure the destination directory exists
mkdir -p "$dest"

# Copy matching .faa files 
for g in "${genomes[@]}"; do
    file="${src}/${g}_scaffolds.proteins.faa"
    if [[ -f "$file" ]]; then
        echo "📂 Copying $file → $dest/"
        cp "$file" "$dest/"
    else
        echo "File not found: $file"
    fi
done

echo "Copy complete!"

####
# the program phyling takes a while, in your best interest to run as a batch script
#### 

#!/bin/bash
#SBATCH --job-name=4gvog_phyling_fungi_aln
#SBATCH --partition=icemhh
#SBATCH --account=icemhh

## 3 day max run time for public partitions, except 4 hour max runtime for the sandbox partition
#SBATCH --time=05-00:00:00 ## time format is DD-HH:MM:SS
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=600G ## max amount of memory per node you require

#SBATCH --error=job%A.err ## %A - filled with jobid
#SBATCH --output=job%A.out ## %A - filled with jobid
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE,TIME_LIMIT_80
#SBATCH --mail-user=swhitner@hawaii.edu

ml lang/Anaconda3/2024.02-1
source activate BUSCO

# directories 
base="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal"
in_dir="$base/fun_proteins"
aln_out="$base/fun_phyling_aln"
tree_out="$base/4gvog_fungi_tree"

mkdir -p "$aln_out" "$tree_out" 

cd "$base"

# 1. Align
echo "Running Phyling alignment..."
phyling align \
  -I /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal/fun_proteins \
  -o /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal/fun_phyling_aln \
  -m fungi_odb10 \
  -t 20 \
  -v

# 2. Build tree 
echo "Building tree..."
phyling tree \
  -I /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal/fun_phyling_aln \
  -o /home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/frank_sugg/test5/my_seqs/fungal/4gvog_fungi_tree \
  -M iqtree \
  -c \
  -p \
  -t 20

echo "Phyling fungal pipeline complete! - fungal tree is ready :) "
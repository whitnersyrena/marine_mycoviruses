########################################################################
#
# This is PART 1 of the contig building plots. 
# this script parses the AUGUSTUS and prodial gff and domtbl files so I can build contig plots in R  
# Syrena Whitner 
# December 2025
#
########################################################################
# for Augustus .gff files files 
########################################################################
#!/bin/bash

# directories 
input_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/syrena_seqs/augustus_output"
output_dir="${input_dir}/parsed_gff"
mkdir -p "$output_dir"

for gff in "${input_dir}"/*.augustus.gff3; do
  base=$(basename "$gff" .augustus.gff3)
  out_csv="${output_dir}/${base}_augustus_features.csv"
  echo "Parsing $base"
  awk -v OFS="," '
  BEGIN { print "contig,feature,start,end,strand,parent" }
  $0 !~ /^#/ {
    contig=$1; feature=$3; start=$4; end=$5; strand=$7; attrs=$9
    split(attrs,a,/[;=]/)
    parent="."
    for(i=1;i<=length(a);i++){
      if(a[i]=="Parent")parent=a[i+1]
      if(a[i]=="ID"&&parent==".")parent=a[i+1]
    }
    if(feature=="gene"||feature=="CDS"||feature=="intron"||feature=="stop_codon"){
      if(feature=="gene"&&prev_contig==contig&&start>prev_end+1){
        print contig,"intergenic",prev_end+1,start-1,".","." 
      }
      print contig,feature,start,end,strand,parent
      if(end>prev_end||prev_contig!=contig){prev_end=end;prev_contig=contig}
    }
  }' "$gff" > "$out_csv"
done

echo "done"

########################################################################
# for prodigal files 
########################################################################
#!/bin/bash

# directories 
input_dir="/home/swhitner/cmaiki_koastore/swhitner/viruses/Sept2025/genomad_output_subset/FINAL/syrena_seqs/hmm_output_common"
output_dir="${input_dir}/parsed_domtbls"
mkdir -p "$output_dir"

# ------------------------------------------------------
# Step 1: Parse each .domtbl into individual CSVs
# ------------------------------------------------------
for f in "${input_dir}"/*.domtbl; do
  base=$(basename "$f" .domtbl)
  out_csv="${output_dir}/${base}_gvog_domains_genomic.csv"
  echo "Parsing $base"
  gawk -v OFS="," '
  BEGIN {
    print "contig,orf_id,feature,alifrom,alito,orf_start,orf_end,strand,dom_genomic_start,dom_genomic_end,evalue,score"
  }
  $0 !~ /^#/ {
    tgt=$1; feat=$4; e=$7; sc=$8; af=$18; at=$19
    if (match($0, /# ([0-9]+) # ([0-9]+) # (-?1) # ID=([^;]+)/, m)) {
      os=m[1]+0; oe=m[2]+0; sd=m[3]; oid=m[4]
      contig=tgt; sub(/_[0-9]+$/, "", contig)
      if (sd==1)  { dgs=os+(af-1)*3; dge=os+at*3-1 }
      else if (sd==-1){ dge=oe-(af-1)*3; dgs=oe-at*3+1 }
      else { dgs=dge="." }
      print contig,oid,feat,af,at,os,oe,(sd==1?"+":"-"),dgs,dge,e,sc
    }
  }' "$f" > "$out_csv"
done

echo "Parsing complete."

# ------------------------------------------------------
# Step 2: Merge all *_gvog_domains_genomic.csv by parent prefix
# ------------------------------------------------------
echo "Merging per-genome CSVs..."

cd "$output_dir" || exit 1

# Collect all unique parent prefixes (everything before "_GVOG")
for prefix in $(ls *_gvog_domains_genomic.csv | sed 's/_GVOG.*//' | sort -u); do
  merged="${prefix}_merged_gvog_domains.csv"
  echo "  → Combining files for $prefix ..."
  
  # Grab header from first file, then append all without headers
  head -n 1 "${prefix}"*_gvog_domains_genomic.csv > "$merged"
  tail -q -n +2 "${prefix}"*_gvog_domains_genomic.csv >> "$merged"
done

echo "All merges done."

####### 
# at this point, you should have two different .csv files! one from augustus and one from prodigal
# these will be the input into the R script titaled "contig_plots.R" in addition to the original protein files for each isolate 

# ---- Paths (relative to /work) ----
CACHE=./vep/cache
DATA=./vep/data
OUT=./vep/out
GVCF_DIR=./output1/variant_calling/haplotypecaller/joint_variant_calling
PLUGINS_DIR=./vep/plugins


# ---- Plugin data files (put them in ./vep/data) ----
REVEL=$DATA/new_tabbed_revel_grch38.tsv.gz
CLINPRED=$DATA/ClinPred_hg38_sorted_tabbed.tsv.gz
LOFTOOL=$DATA/LoFtool_scores.txt
AMISS=$DATA/AlphaMissense_hg38.tsv.gz
DBNSFP=$DATA/dbNSFP4.9c_GRCh38.tsv.bgz

# index helper
ensure_tbi () {
  local f="$1"
  if [[ -f "$f" && ! -f "$f.tbi" ]]; then
    echo "Indexing $f with tabix..."
    tabix -s 1 -b 2 -e 2 "$f"
  fi
}

# build plugin args conditionally (only add if file exists)
PLUGIN_ARGS=()
if [[ -f "$REVEL" ]];   then PLUGIN_ARGS+=(--plugin "REVEL,$REVEL"); fi
if [[ -f "$CLINPRED" ]];then PLUGIN_ARGS+=(--plugin "ClinPred,$CLINPRED"); fi
if [[ -f "$LOFTOOL" ]]; then PLUGIN_ARGS+=(--plugin "LoFtool,$LOFTOOL"); fi
if [[ -f "$AMISS" ]];   then PLUGIN_ARGS+=(--plugin "AlphaMissense,file=$AMISS"); fi
if [[ -f "$DBNSFP" ]];  then PLUGIN_ARGS+=(--plugin "dbNSFP,$DBNSFP,SIFT_score,Polyphen2_HDIV_score,MutationTaster_score"); fi


# common VEP options
VEP_COMMON=(
  --offline --cache --dir_cache "$CACHE" --assembly GRCh38
  --dir_plugins "$PLUGINS_DIR"
  --species homo_sapiens
  --vcf --compress_output bgzip
  --fork 4
)

# ---- annotate all per-sample gVCFs ----
shopt -s nullglob
GVCFS=("$GVCF_DIR"/*.vcf.gz)
if (( ${#GVCFS[@]} > 0 )); then
  echo "Annotating ${#GVCFS[@]} per-sample gVCFs..."
  for V in "${GVCFS[@]}"; do
    BN=$(basename "$V" .vcf.gz)
    OUTVCF="$OUT/${BN}.vep.vcf.gz"
    echo "  -> $OUTVCF"
    vep \
      "${VEP_COMMON[@]}" \
      --input_file "$V" \
      --output_file "$OUTVCF" \
      "${PLUGIN_ARGS[@]}"
  done
else
  echo "No per-sample gVCFs found in $GVCF_DIR"
fi

echo "Done. Outputs in $OUT"
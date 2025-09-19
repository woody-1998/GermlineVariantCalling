# ---- Paths (relative to /work) ----
GV_DIR=./gvanno
GV_PY=$GV_DIR/gvanno.py
INPUT=./output1/variant_calling/haplotypecaller/joint_variant_calling
OUT=$GV_DIR/out

mkdir -p "$OUT"


# ---- use gvanno to annotate the file----
shopt -s nullglob
VCFS=("$INPUT"/*.vcf.gz)
if (( ${#VCFS[@]} > 0 )); then
  echo "annotating ${#VCFS[@]} per-sample VCFs..."
  for V in "${VCFS[@]}"; do
    BN=$(basename "$V" .vcf.gz)
    echo "now the sample is ${BN}"
    python3 $GV_PY \
      --query_vcf "$V" \
      --gvanno_dir "$GV_DIR" \
      --output_dir "$OUT" \
      --genome_assembly grch38 \
      --sample_id "${BN}" \
      --container docker \
      --oncogenicity_annotation \
      --vep_lof_prediction \
      --vep_gencode_basic \
      --force_overwrite

  done
else
  echo "No per-sample gVCFs found in $VCF_DIR"
fi

echo "Done. Outputs in $OUT"
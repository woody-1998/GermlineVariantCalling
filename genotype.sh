# ---- Paths (relative to /work) ----
REF=./ref/Homo_sapiens/GATK/GRCh38/Sequence/WholeGenomeFasta/Homo_sapiens_assembly38.fasta
GVCF_DIR=./output1/variant_calling/haplotypecaller
OUT=./output1/variant_calling/genotyped


# ---- genotype per-sample gVCFs to VCF files----
shopt -s nullglob
GVCFS=("$GVCF_DIR"/*/*.g.vcf.gz)
if (( ${#GVCFS[@]} > 0 )); then
  echo "genotyping ${#GVCFS[@]} per-sample gVCFs..."
  for V in "${GVCFS[@]}"; do
    BN=$(basename "$V" .g.vcf.gz)
    OUTVCF="$OUT/${BN}.vcf.gz"
    echo "  -> $OUTVCF"
    gatk GenotypeGVCFs \
      -R "$REF" \
      -V "$V" \
      -O "$OUTVCF" \

  done
else
  echo "No per-sample gVCFs found in $GVCF_DIR"
fi

echo "Done. Outputs in $OUT"
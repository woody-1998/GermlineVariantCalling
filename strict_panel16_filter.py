#!/usr/bin/env python3
"""
strict_panel16_filter.py

Strictly filter a joint gvanno TSV to mirror a Table-2-style pathogenic set,
restricted to the 16 hereditary cancer genes. Produces a one-row-per-gene
summary including: Gene, Mutation type, rsID, Impact, Cancer history, Exon,
plus HGVS c/p and n_probands.

Usage:
  python3 strict_panel16_filter.py \
    --in_tsv /path/joint_germline_recalibrated_gvanno_grch38.pass.tsv \
    --out_tsv /path/panel16_STRICT_one_per_gene_with_details.tsv \
    [--af_cutoff 1e-4] \
    [--min_stars 2] \
    [--keep_noncanonical] \
    [--exon_mode raw|number]    # raw = e.g. "16/16"; number = "16" (default)

Notes:
- ClinVar filter: accept tokens "Pathogenic" or "Likely_pathogenic"; reject
  benign/conflicting/risk/association/protective/drug_response/affects/uncertain.
- Review status: stars >= min_stars (default 2).
- Rarity: gnomAD AF < af_cutoff (default 1e-4) or missing.
- Consequence: keep LoF (frameshift, stop_gained, splice_acceptor/donor, start/stop_lost).
  Special-case: keep MUTYH p.Gly396Asp (c.1187G>A, rs36053993).
- Canonical: default ON; use --keep_noncanonical to disable.
"""

import argparse, os, re
import pandas as pd

PANEL16 = {"APC","ATM","BRCA1","BRCA2","CDH1","CHEK2","EPCAM","MLH1","MSH2","MSH6",
           "PMS2","PTEN","PALB2","STK11","TP53","MUTYH"}

LOF_TERMS = ("frameshift_variant","stop_gained","splice_acceptor_variant",
             "splice_donor_variant","start_lost","stop_lost")

def find_col(df: pd.DataFrame, candidates):
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None

def parse_args():
    ap = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("--in_tsv", required=True, help="Joint gvanno TSV (.tsv or .tsv.gz)")
    ap.add_argument("--out_tsv", required=True, help="Output TSV path")
    ap.add_argument("--af_cutoff", type=float, default=1e-4, help="Max gnomAD AF to keep (exclusive)")
    ap.add_argument("--min_stars", type=int, default=2, help="Minimum ClinVar review stars to keep")
    ap.add_argument("--keep_noncanonical", action="store_true", help="Do not restrict to canonical transcripts")
    ap.add_argument("--exon_mode", choices=["raw","number"], default="number",
                    help="Output Exon as raw (e.g., '16/16') or numeric-only (e.g., '16')")
    return ap.parse_args()

def main():
    args = parse_args()

    df = pd.read_csv(args.in_tsv, sep="\t", dtype=str)

    # Map columns (case-insensitive + common aliases)
    SYMBOL = find_col(df, ["SYMBOL","Gene","gene_symbol"]) or "SYMBOL"
    HGVSc  = find_col(df, ["HGVSc","HGVS_c","HGVSc_change"]) or "HGVSc"
    HGVSp  = find_col(df, ["HGVSp","HGVS_p","HGVSp_change"]) or "HGVSp"
    CLNSIG = find_col(df, ["CLINVAR_CLNSIG","clinvar_clnsig","clinvar_significance"]) or "CLINVAR_CLNSIG"
    STARS  = find_col(df, ["CLINVAR_REVIEW_STATUS_STARS","clinvar_review_status_stars"]) or "CLINVAR_REVIEW_STATUS_STARS"
    GNOMAD = find_col(df, ["gnomADg_AF","gnomAD_AF","gnomadg_af"]) or "gnomADg_AF"
    IMPACT = find_col(df, ["IMPACT","impact"]) or "IMPACT"
    CONSEQ = find_col(df, ["CONSEQUENCE","Consequence","consequence"]) or "CONSEQUENCE"
    EXON   = find_col(df, ["EXON","Exon","exon"]) or "EXON"
    SAMPLE = find_col(df, ["SAMPLE_ID","sample_id"]) or "SAMPLE_ID"
    CANON  = find_col(df, ["CANONICAL","canonical"])
    TRAITS = find_col(df, ["CLINVAR_TRAITS_ALL","clinvar_traits_all"])
    EXIST  = find_col(df, ["Existing_variation","ID","dbsnp_rsid"])

    # Restrict to 16-gene panel
    x = df[df[SYMBOL].isin(PANEL16)].copy()

    # ClinVar P/LP only (reject uncertain/benign/etc.)
    sig = x[CLNSIG].astype(str).str.lower()
    accept = sig.str.contains(r'(?:^|[|,;:_\s])(likely_pathogenic|pathogenic)(?:$|[|,;:_\s])', regex=True)
    reject = sig.str.contains(r'benign|conflicting|risk_factor|association|protective|drug_response|affects|uncertain', regex=True)
    x = x[accept & ~reject].copy()

    # Review stars >= min_stars
    def stars_ok(v):
        if pd.isna(v): return False
        m = re.search(r'(\d+)', str(v))
        return int(m.group(1)) >= args.min_stars if m else False
    x = x[x[STARS].apply(stars_ok)].copy()

    # AF < cutoff or missing
    def af_ok(v):
        if pd.isna(v) or v==".": return True
        try: return float(v) < args.af_cutoff
        except: return True
    x = x[x[GNOMAD].apply(af_ok)].copy()

    # Canonical transcript only (unless user disables)
    if CANON and not args.keep_noncanonical:
        x = x[x[CANON].astype(str).str.upper()=="YES"].copy()

    # Keep LoF + allow MUTYH p.Gly396Asp
    is_lof = x[CONSEQ].astype(str).apply(lambda s: any(t in s for t in LOF_TERMS))
    is_mutyh_known = (x[SYMBOL]=="MUTYH") & (
        x[HGVSp].astype(str).str.contains("Gly396Asp") |
        x[HGVSc].astype(str).str.contains(r"c\.1187G>A") |
        (x[EXIST].astype(str).str.contains(r"\brs36053993\b", case=False, regex=True) if EXIST else False)
    )
    x = x[is_lof | is_mutyh_known].copy()

    # Derive Exon + rsID
    if args.exon_mode == "number":
        def exon_num(s):
            if pd.isna(s): return ""
            m = re.search(r'(\d+)(?:/\d+)?', str(s))
            return m.group(1) if m else ""
        x["Exon"] = x[EXON].apply(exon_num)
    else:
        x["Exon"] = x[EXON].fillna("")

    def pick_rsid(v):
        if v is None: return "–"
        m = re.search(r'(rs\d+)', str(v), re.I)
        return m.group(1) if m else "–"
    x["rsID"] = x[EXIST].apply(pick_rsid) if EXIST else "–"

    # Collapse to locus-level & count carriers
    grp_cols = ["CHROM","POS","REF","ALT",SYMBOL]
    agg = (x.groupby(grp_cols, dropna=False)
             .agg({
                HGVSc:'first', HGVSp:'first', CONSEQ:'first', IMPACT:'first',
                "Exon":'first', "rsID":'first', (TRAITS if TRAITS else CONSEQ):'first',
                SAMPLE: pd.Series.nunique
             })
             .reset_index()
             .rename(columns={SAMPLE:"n_probands", SYMBOL:"Gene",
                              HGVSc:"HGVS_DNA", HGVSp:"HGVS_protein",
                              CONSEQ:"Mutation type", IMPACT:"Impact"}))

    if TRAITS and TRAITS in agg.columns:
        agg = agg.rename(columns={TRAITS:"Cancer history"})
    else:
        agg["Cancer history"] = ""

    # One row per gene: prioritize frameshift/stop_gained > splice > start/stop_lost, then higher n_probands
    def lof_priority(mut_type: str) -> int:
        s = str(mut_type or "")
        if "frameshift_variant" in s or "stop_gained" in s: return 3
        if "splice_acceptor_variant" in s or "splice_donor_variant" in s: return 2
        if "start_lost" in s or "stop_lost" in s: return 1
        return 0

    agg["rank"] = agg["Mutation type"].apply(lof_priority)
    one = (agg.sort_values(["Gene","rank","n_probands"], ascending=[True, False, False])
              .groupby("Gene", as_index=False).first())[slice(0, 8)]

    # Final columns
    final_cols = ["Gene","Mutation type","rsID","Impact","Cancer history","Exon",
                  "HGVS_DNA","HGVS_protein","n_probands"]
    one[final_cols].to_csv(args.out_tsv, sep="\t", index=False)
    print(f"Wrote {args.out_tsv} with {len(one)} genes. Exon mode: {args.exon_mode}")

if __name__ == "__main__":
    main()

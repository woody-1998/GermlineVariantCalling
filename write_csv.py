import os
import csv
import glob



fa_list = glob.glob("./germline/*.fastq.gz")
out_csv = "testsheet.csv"
samples = {}

print(fa_list)
print(len(fa_list))

for fa in fa_list:
    base = os.path.basename(fa)
    run = base.split("_")[0] # ERR number prefix
    if run not in samples:
        samples[run] = {}
    if "_1" in base:
        samples[run]["f1"] = fa
    elif "_2" in base:
        samples[run]["f2"] = fa

print(samples)


with open(out_csv, "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["patient", "sample", "lane", "fastq_1", "fastq_2"])
    for i, (run, d) in enumerate(samples.items(), 1):
        if i >= 2 :
            break
        if "f1" in d and "f2" in d:
            patient = run  # or use f"ID{i}"
            sample = f"S{i}"
            lane = "L001"
            writer.writerow([patient, sample, lane, d["f1"], d["f2"]])

print(f"Wrote samplesheet: {out_csv}")

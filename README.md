This is a germline variant calling workflow trying to repeat the [Study of Germline Varaints in Nothern Brazil](https://link.springer.com/content/pdf/10.1186/s12885-021-08089-9.pdf)

1. Conda environment:
   conda config --add channels bioconda
   conda create -n nextflow python pandas matplotlib bcftools nextflow gatk4
   conda activate nextflow
   brew install awscli

2. Nextflow/nf-core/sarek
   * download genome reference documents under the directory ""./ref" :
     aws s3 cp --no-sign-request --recursive s3://ngi-igenomes/igenomes/Homo_sapiens/GATK/GRCh38/ ./ref/Homo_sapiens/GATK/GRCh38/
     
   * run the write_csv.py get the input csv file
     
   * run the piprline:
     nextflow run nf-core/sarek  -profile docker \
       --input samplesheet1.csv \
       --outdir ./output1  \
       --tools haplotypecaller \
       --joint_germline true  \
       -c resourses.config  \
       --igenomes_base <absolute_path_to_ref>  \
       -resume
     
   * genotyping GVCF files to VCF files
     bash genotype.sh (check the path)

3. Vep Annotation
   * pull the docker
     docker pull ensembl-vep:release_115.0
     
   * install and download all plugins following the [official manual](https://www.ensembl.org/info/docs/tools/vep/script/vep_plugins.html#revel)
  
   * docker run --rm -it -v $PWD:/work ensemblorg/ensembl-vep:release_115.0 \
     cd /work \
     bash annotation.sh

4. Gvanno Annotation
   * pull the docker
     docker pull sigven/gvanno:1.7.0

   * docker run --rm -it -v $PWD:/work sigven/gvanno:1.7.0 \
     cd \work
     bash gvanno.sh
     
  

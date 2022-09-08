#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

// container -> /lustre/scratch118/humgen/resources/ensembl/vep/singularity_containers/vep_104.0.sif

params.dir_cache = "/lustre/scratch118/humgen/resources/ensembl/vep/GRCh38/vep_data"
params.fasta = "/lustre/scratch118/humgen/resources/ensembl/vep/GRCh38/vep_data/homo_sapiens/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
params.plugins = "/lustre/scratch118/humgen/resources/ensembl/vep/GRCh38/Plugins"
params.vcf_file = "./testdata.vcf"
params.cad_path = "/lustre/scratch118/humgen/resources/cadd_scores/20201027-GRCh38_v1.6"
params.spliceAI = "/lustre/scratch118/humgen/resources/SpliceAI_data_files" 
params.fork = 4

process vep {
  input: 
    val x
  output:
    stdout
  script:
    """
vep \
--cache \
--dir_cache $params.dir_cache \
--fasta $params.fasta \
--offline \
--format vcf \
--dir_plugins $params.plugins \
-i $params.vcf_file \
--plugin SpliceRegion,Extended \
--plugin GeneSplicer,$params.plugins/GeneSplicer/bin/linux/genesplicer,$params.plugins/GeneSplicer/human \
--plugin UTRannotator,$params.plugins/uORF_starts_ends_GRCh38_PUBLIC.txt \
--plugin CADD,$params.cad_path/whole_genome_SNVs.tsv.gz,$params.cad_path/gnomad.genomes.r3.0.indel.tsv.gz \
--plugin dbNSFP,$params.plugins/dbNSFP4.1a.gz,Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
--plugin DisGeNET,file=$params.plugins/all_variant_disease_pmid_associations_final.tsv.gz \
--fork $params.fork \
--everything \
--plugin Phenotypes,file=$params.plugins/phenotypes.gff.gz,include_types=Variation \
--plugin Conservation,$params.plugins/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw \
--plugin LoF,loftee_path:$params.plugins,human_ancestor_fa:$params.plugins/GRCh38_human_ancestor.fa.gz,conservation_file:$params.plugins/loftee.sql,gerp_bigwig:$params.plugins/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw,debug:1 \
--plugin REVEL,$params.plugins/grch38_tabbed_revel.tsv.gz \
--plugin SpliceAI,snv=$params.spliceAI/spliceai_scores.raw.snv.hg38.vcf.gz,indel=$params.spliceAI/spliceai_scores.raw.indel.hg38.vcf.gz \
--vcf \
-o /opt/vcf/$OUTPUT_VCF \
--compress_output bgzip \
--allele_number \
--verbose 
    """
}

workflow {
  Channel.of('test') | vep | view
}
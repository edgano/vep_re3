#!/usr/bin/env nextflow
nextflow.enable.dsl=2 

// container -> /lustre/scratch118/humgen/resources/ensembl/vep/singularity_containers/vep_104.0.sif

// Directories to BIND
dir_cache = "/lustre/scratch118/humgen/resources/ensembl/vep/GRCh38/vep_data"
plugins = "/lustre/scratch118/humgen/resources/ensembl/vep/GRCh38/Plugins"
exomes = "/lustre/scratch118/humgen/resources/gnomAD/release-2.1.1/exomes"
cad_path = "/lustre/scratch118/humgen/resources/cadd_scores/20201027-GRCh38_v1.6"
spliceAI = "/lustre/scratch118/humgen/resources/SpliceAI_data_files"

bindSite = "/opt/vep/.vep"
bindSiteCache = "${bindSite}"
bindSitePlugins = "${bindSite}/Plugins"

// variables
params.fasta = "${bindSite}/homo_sapiens/104_GRCh38/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
params.vcf_file = "/nfs/users/nfs_e/en6/vep_re3/testdata.vcf"
params.snv = "/lustre/scratch118/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.snv.hg38.vcf.gz"
params.indel= "/lustre/scratch118/humgen/resources/SpliceAI_data_files/spliceai_scores.raw.indel.hg38.vcf.gz"
params.compress_output="bgzip"
UTRannotatorInput = "${bindSitePlugins}/uORF_starts_ends_GRCh38_PUBLIC.txt"
cadSNVs="/lustre/scratch118/humgen/resources/cadd_scores/20201027-GRCh38_v1.6/whole_genome_SNVs.tsv.gz"
cadIndel="/lustre/scratch118/humgen/resources/cadd_scores/20201027-GRCh38_v1.6/gnomad.genomes.r3.0.indel.tsv.gz"
dbNSFPinput="${bindSitePlugins}/dbNSFP4.1a.gz"
DisGeNETinput="${bindSitePlugins}/all_variant_disease_pmid_associations_final.tsv.gz"
phenotypesFile="${bindSitePlugins}/phenotypes.gff.gz"
conservationScore="${bindSitePlugins}/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw"
humanAncesforFasta="${bindSitePlugins}/GRCh38_human_ancestor.fa.gz"
conservationFile="${bindSitePlugins}/loftee.sql"
gerpFile="${bindSitePlugins}/90_mammals.gerp_conservation_scores.homo_sapiens.GRCh38.bw"
inputRebel= "${bindSitePlugins}/grch38_tabbed_revel.tsv.gz"
params.fork = 4

process vep {
    containerOptions "--bind ${dir_cache}:${bindSiteCache} \
                    --bind ${plugins}:${bindSitePlugins} \
                    --bind ${exomes} \
                    --bind ${cad_path} \
                    --bind ${spliceAI} "
  input:
    file vcf

  output:
    path "variant_effect_output.*"

  script:
    """
    vep \
    --cache \
    --dir_cache ${bindSiteCache}/ \
    --fasta $params.fasta \
    --offline \
    --format vcf \
    --dir_plugins ${bindSitePlugins} \
    -i $vcf \
    --plugin SpliceRegion,Extended \
    --plugin GeneSplicer,${bindSitePlugins}/GeneSplicer/bin/linux/genesplicer,${bindSitePlugins}/GeneSplicer/human \
    --plugin UTRannotator, ${UTRannotatorInput}\
    --plugin CADD,${cadSNVs},${cadIndel}\
    --plugin dbNSFP,${dbNSFPinput},Ensembl_transcriptid,Uniprot_acc,VEP_canonical,LRT_pred,SIFT_pred,MutationTaster_pred,Polyphen2_HDIV_pred,Polyphen2_HVAR_pred \
    --plugin DisGeNET,file=${DisGeNETinput} \
    --fork ${params.fork} \
    --everything \
    --plugin Phenotypes,file=${phenotypesFile},include_types=Variation \
    --plugin Conservation,${conservationScore} \
    --plugin LoF,loftee_path:${bindSitePlugins},human_ancestor_fa:${humanAncesforFasta},conservation_file:${conservationFile},gerp_bigwig:${gerpFile},debug:1 \
    --plugin REVEL,${inputRebel} \
    --plugin SpliceAI,snv=${params.snv},indel=${params.indel} \
    --vcf \
    --compress_output ${params.compress_output} \
    --allele_number \
    --verbose
    """
}

workflow {
    vcf_ch = Channel.fromPath("$params.vcf_file")

    vep(vcf_ch)
}
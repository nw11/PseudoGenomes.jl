include(Pkg.dir("PseudoGenomes", "src","pseudogenome-snp.jl"))
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-snp.jl"))

mm10path =joinpath(homedir(),"Dropbox/Data/GENOME_DATA/mm10full/mm10.fa")
mm10_masked_fasta = joinpath( tempdir(), "mm10.castmaskdeletion.fasta")

SVs_filename = Pkg.dir( "PseudoGenomes","testdata","28strains.REL-1410-SV.sdp.SAMPLE.tab.gz")
mask_fasta_file_with_strain_SV_deletions(mm10path,mm10_masked_fasta,SVs_filename,"CAST_EiJ")

mm10_masked_and_snp_subst_fasta = joinpath( tempdir(), "mm10.castmaskdeletion.snpsubst.fasta")
substitue_snps_in_fasta_file(mm10_masked_fasta,vcf_path,mm10_masked_and_snp_subst_fasta)

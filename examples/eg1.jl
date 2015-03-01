using FactCheck
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-snp.jl"))
mm10path =joinpath(homedir(),"Dropbox/Data/GENOME_DATA/mm10full/mm10.fa")
vcf_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata", "vcf-file1.vcf")
outpath = joinpath( Pkg.dir(),"PseudoGenomes","testdata","out.fasta")
substitue_snps_in_fasta_file(mm10path,vcf_path,outpath)

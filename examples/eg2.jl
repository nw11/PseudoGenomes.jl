# Mask deletions
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-deletions.jl"))
mm10path =joinpath(homedir(),"Dropbox/Data/GENOME_DATA/mm10full/mm10.fa")
outpath = joinpath( tempdir(), "mm10.castmaskdeletion.fasta")
SVs_filename = Pkg.dir( "PseudoGenomes","testdata","28strains.REL-1410-SV.sdp.SAMPLE.tab.gz")
mask_fasta_file_with_strain_SV_deletions(mm10path,outpath,SVs_filename,"CAST_EiJ")


using FactCheck
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-deletions.jl"))

ref_fasta = joinpath( Pkg.dir(),"PseudoGenomes","testdata","test.fasta")

# remember to change the coordinate strings in the strain fields if you want to play around with length of deletion masked.
SVs_filename = Pkg.dir("PseudoGenomes","testdata","28strains.REL-1410-SV.sdp.SAMPLE.MOD.tab")
outpath = joinpath( Pkg.dir(),"PseudoGenomes","testdata","out.maskdeletions.fasta")

# -- test
SV_df = read_SV_file(SVs_filename, gzip=false)

# -- test
#"PWK_PhJ"
strain_SVs_df = get_strain_SV_type(SV_df, "129P2_OlaHsd",Regex("DEL"))

# -- test
mask_fasta_file_with_strain_SV_deletions(ref_fasta,outpath,strain_SVs_df)

# -- test
#mask_region!(sequence,SV_deletions_df)

#mask_fasta_file_with_strain_SV_deletions(mm10path,outpath,strain_SVs_df)
#mask_fasta_file_with_strain_SV_deletions(mm10path,outpath,SVs_filename,"CAST_EiJ")

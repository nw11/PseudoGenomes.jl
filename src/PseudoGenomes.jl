module PseudoGenomes

# package code goes here

export substitute_snps_in_fasta_file
export check_variant_alleles_in_fasta_file
export read_snp_positions_from_columns
export read_vcf_columns

include(Pkg.dir("PseudoGenomes","src","pseudogenome-snp.jl"))

export mask_fasta_file_with_strain_SV_deletions
export get_strain_SVs

include(Pkg.dir("PseudoGenomes","src","pseudogenome-deletions.jl"))


export make_diploid
include(Pkg.dir("PseudoGenomes","src","pseudogenome-deletions.jl"))
end # module

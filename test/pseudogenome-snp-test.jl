using FactCheck
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-snp.jl"))

ref_fasta = joinpath( Pkg.dir(),"PseudoGenomes","testdata","test.fasta")
vcf_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata", "vcf-file2.vcf")
outpath = joinpath( Pkg.dir(),"PseudoGenomes","testdata","out.fasta")

# - test substitute_
expected_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata","expected.fasta")
substitute_snps_in_fasta_file(ref_fasta,vcf_path,outpath)

expected =  readall(expected_path)
got = readall(outpath)

facts("fasta-check") do
    @fact got => expected
end

# - test check variant alleles
check_variant_alleles_in_fasta_file(outpath,vcf_path)


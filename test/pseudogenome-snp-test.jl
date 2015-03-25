using FactCheck
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-snp.jl"))

ref_fasta = joinpath( Pkg.dir(),"PseudoGenomes","testdata","test.fasta")
vcf_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata", "vcf-file2.vcf")
vcf_noheader_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata", "vcf-file2.noheader.vcf")
outpath = joinpath( Pkg.dir(),"PseudoGenomes","testdata","out.fasta")
cpgbedfile1 = joinpath( Pkg.dir(),"PseudoGenomes","testdata","cpgs1.bedgraph")
#cpgbedfile2 = joinpath( Pkg.dir(),"PseudoGenomes","testdata","cpgs2.bedgraph")
# - test read_vcf_columns
read_snp_positions_from_columns(vcf_path, gzip=false )

# - test substitute_
expected_path = joinpath( Pkg.dir(),"PseudoGenomes","testdata","expected.fasta")
substitute_snps_in_fasta_file(ref_fasta,vcf_path,outpath)

expected =  readall(expected_path)
got = readall(outpath)

facts("fasta-check") do
    @fact got => expected
end


# -  test cpg
a = ['A','C','C','C','G','T','G','C','G','A']
cgpos = find_cpgs(a)
# default zero based
expected_cgpos = [3,7]
facts("cpg-check") do
    @fact cgpos => expected_cgpos
end


# - make cpg bedfile
#make_cpg_bedfile_type1(outpath,cpgbedfile1)
make_cg_bedfile(outpath,cpgbedfile1)


# - test check variant alleles
# - switch logging to file and check file is correct
check_variant_alleles_in_fasta_file(outpath,vcf_path)

# -
# -
check_variant_alleles_in_fasta_file(outpath,vcf_noheader_path)

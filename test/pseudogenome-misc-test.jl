using FactCheck
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-misc.jl"))

fastafile1 = Pkg.dir("PseudoGenomes","testdata","test.fasta")
fastafile2 = Pkg.dir("PseudoGenomes","testdata","expected.fasta")
outputfasta = Pkg.dir("PseudoGenomes","testdata","diploid.fasta")
id1="abc"
id2="def"

make_diploid(fastafile1, id1, fastafile2, id2, outputfasta)

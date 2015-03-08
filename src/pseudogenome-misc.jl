using FastaIO
using Lumberjack
function ens2ucsc(seq_id::String)
    if seq_id == "MT"
        return "chrM"
    else
        return "chr"*seq_id
    end
end

function ucsc2ens(seq_id::String)
   if seq_id == "chrM"
        return "MT"
    else
        return "chr"*seq_id
    end
end

# yes this could be done via sed and cat
# but here it is.
function make_diploid(fastafile1, id1, fastafile2, id2, output_fasta)
    #run through each assign id to description
    fw = FastaWriter(output_fasta)
    fr1 = FastaReader{Vector{Char}}(fastafile1)
    Lumberjack.info("Start reading fasta file $fastafile1")
    for (desc, seq) in fr1
        desc = desc * "_$id1"
        writeentry(fw, desc, seq)
    end
    close(fr1)

    fr2 = FastaReader{Vector{Char}}(fastafile2)
    Lumberjack.info("Start reading fasta file $fastafile2")
    for (desc, seq) in fr2
        desc = desc * "_$id2"
        writeentry(fw, desc, seq)
    end
    close(fr2)
    close(fw)
    Lumberjack.info("Done writing fasta file $output_fasta")
end


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

function make_diploid(fastafile1, id1, fastafile2, id2, output_fasta)
    #run through each assign id to description
end

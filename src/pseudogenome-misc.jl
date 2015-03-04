
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

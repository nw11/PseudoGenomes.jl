using FastaIO
using DataFrames
using Lumberjack
using GZip

include(Pkg.dir("PseudoGenomes", "src","pseudogenome-misc.jl"))
#
# This is a file formatted in the way you get it
# from WTSI see here: ftp://ftp-mouse.sanger.ac.uk/current_svs/
#
function read_SV_file(filename::ASCIIString;gzip=true)
    num_lines=0
    heading = Regex("^#CHROM")
    line_itr = nothing
    if gzip
       line_itr = eachline(GZip.open(filename) )
    else
       line_itr = eachline(open(filename) )
    end

    # check for comment lines and header
    df = DataFrame()
    header=[]
    for line in line_itr
       line=chomp(line)
       if ismatch(heading, line)
           println("HEADING: $line")
           header = split(line,'\t')
           species_names=header[5:end]
           println("The species names are $species_names")
           break
       end
       # anything above the line heading is a comment
       println("COMMENT: $line")
    end

    first_data_line = first(line_itr)
    first_data_line=chomp(first_data_line)
    first_data_fields = split(first_data_line,'\t')
    df[:seq_id]     = first_data_fields[1]
    df[:start]   = first_data_fields[2]
    df[:stop]    = first_data_fields[3]
    df[:formats] = first_data_fields[4]
    for field_idx = 5:length(header)
        df[symbol(header[field_idx])] = first_data_fields[field_idx]
    end

    println("dataframe after first dataline: $df")
    df_num_cols = ncol(df)
    for line in line_itr
         line=chomp(line)
         fields = split(line,'\t')
         num_fields = length(fields)
         if num_fields == df_num_cols
             push!(df,fields)
         else
             error("Number of fields parsed inconsistent $fields from $line.\nNumber fields parsed: $num_fields Num Cols Expected $df_num_cols")
         end
         num_lines += 1
    end
    #DELS = SVs[ SVs[:SV_type] .== "DEL", :]
    println( "Total number of lines: $num_lines")

    return df
end


function get_strain_SVs(SVs::DataFrame, strain_id::ASCIIString,sv_type::Regex; seq_id_format="ucsc")

    # get the strain
    strain_column = SVs[symbol(strain_id)]

    chrs=ASCIIString[]
    starts=Int64[]
    stops=Int64[]

    # and parse the first field to get the actual coordinates for the structural vairant
    for row in strain_column
        # The fields are POS::CLASSIFICATION::BREAKPPOINT::REPEAT_TYPE
        # it seems there used to be a field called PEM_pattern that crops up.
        # The fields are specified  in the format field of the tab file for each row (POS:CL:BP:TY)
        # so this prob should be take this into account in case they change the fields in the future.

        fields = split(row, ';')

        if length(fields) == 1
            continue
        elseif length(fields) < 4
            Lumberjack.warn("Less than expected number of fields")
        end

        if !ismatch(sv_type,fields[2])
            continue
        end
        coordinates=fields[1]
        (chr, position) = split(coordinates,':')
        (start,stop)=split(position,'-')
        start_int = parseint(start)
        stop_int  = parseint(stop)
        if seq_id_format == "ucsc"
            push!(chrs, ens2ucsc(chr) )
        else
            push!(chrs, chr)
        end
        push!(starts,start_int)
        push!(stops,stop_int)
    end
    return DataFrame(seq_id=chrs,start=starts,stop=stops)
end

function mask_region!(sequence,SV_deletions_df)
   num_masked=0
   for row=1:nrow(SV_deletions_df)
        startpos = SV_deletions_df[row,:start]
        stoppos  = SV_deletions_df[row,:stop]
        for i=startpos:stoppos
            sequence[i]='N'
            num_masked +=1
        end
    end
    Lumberjack.info("Done masking")
    return num_masked
end

function mask_fasta_file_with_strain_SV_deletions(filename,outfilename,SV_deletions::DataFrame)
    # roll through the fastafile and then roll through the fastafile switch to N.

    fr = FastaReader{Vector{Char}}(filename)
    fw = FastaWriter(outfilename)
    Lumberjack.info("Start masking of fasta file")
    for (desc, seq) in fr
         seqlength = length(seq)
         Lumberjack.info("processing $desc - length: $seqlength")
         seqid_df = SV_deletions[ SV_deletions[:seq_id] .== desc,:]
         if nrow(seqid_df) > 0
             num_masked = mask_region!(seq, seqid_df )
             Lumberjack.info("Masked $num_masked nucleotides")
         else
             Lumberjack.info("No deletion coordinates for $desc")
         end

         writeentry(fw, desc, seq)
    end
   Lumberjack.info("read $(fr.num_parsed) entries")
   close(fw)
end

function mask_fasta_file_with_strain_SV_deletions(fastafile::ASCIIString,outfilename::ASCIIString, SVfile::ASCIIString,strain::ASCIIString)
    SV_df = read_SV_file(SVfile)
    strain_SVs_df = get_strain_SVs(SV_df, strain,Regex("DEL"))
    mask_fasta_file_with_strain_SV_deletions(fastafile,outfilename,strain_SVs_df)
end

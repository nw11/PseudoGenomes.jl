using FastaIO
using DataFrames
using Lumberjack
include(Pkg.dir("PseudoGenomes", "src","pseudogenome-deletions.jl"))

function read_snp_positions_from_columns(filename::String; seq_id_format="ucsc", gzip=false, header=true )
    #locate comments and header
    heading_rgx = Regex("^#CHROM")
    line_itr = nothing
    line_array=Any[]
    seq_id = ASCIIString[]
    positions = Int64[]
    ref = Char[]
    variant = Char[]
    line_num = 0
    df = DataFrame()

    if gzip
       #line_itr = eachline(GZip.open(filename) ) - GZip is broken for some versions of gzip - gives an extraline - see bug report
       line_itr  = eachline(`zcat $filename`)
    else
       line_itr  = eachline(open(filename) )
    end


    # check for comment lines and header
    if header
        for line in line_itr
            line=chomp(line)
            if ismatch(heading_rgx, line)
                Lumberjack.info("DETECT HEADING: $line")
                break
            end
            # anything above the line heading is a comment
           Lumberjack.info("DETECT COMMENT: $line")
        end
    end

    for line in line_itr
        line_num += 1
        # split and choose columns
        line_array =  map(x->strip(x), split(line,'\t') )

        if seq_id_format == "ucsc"
            push!(seq_id, ens2ucsc( line_array[1]) )
        else
           push!(seq_id, line_array[1] )
        end
        push!(positions, parseint( line_array[2]))
        push!(ref,line_array[4][1] )
        push!(variant,line_array[5][1] )
        if (line_num % 1000000) == 0
            Lumberjack.info("read $line_num lines")
        end
    end
    return DataFrame(seq_id = seq_id, position = positions, ref=ref, variant=variant)
end

function read_vcf_columns(filename::String, cols, sequence_id_col=1; gzip="false" )
    sequence_dict = Dict()
    ignore_rgx = Regex("^#")
    for line in eachline(open(filename))
       if ismatch(ignore_rgx, line)
           continue
       end
       # split and choose columns
       line_array =  map(x->strip(x), split(line,'\t')[cols])
       seq_id = line_array[sequence_id_col]

       if ! haskey( sequence_dict, seq_id )
           println("Don't see-",seq_id)
           df = DataFrame()
           for i = 1:length(cols)
               df[i] = Any[]
           end
           sequence_dict[ seq_id ] =df
       end

       # push onto the dataframe
       push!(sequence_dict[   seq_id  ],line_array)
    end
    return sequence_dict
end

# See vcf-file1.vcf for the kind of file you would expect this function to take.
# Essentially a vcf file or one without a heading, but with the mandatory columns present.
# We don't have a vcf parser, and dataframes readtable doesnt allow to choose
# columns meaning we don't get it's type inferenc as it would do the job.
function read_snp_positions_from_columnsOld(filename::String,seq_id_format="ucsc" )
    sequence_dict = Dict()
    ignore_rgx = Regex("^#")
    cols = [1,2,4,5]
    line_array=Any[]
    seq_id = ASCIIString[]
    positions = Int64[]
    ref = Char[]
    variant = Char[]
    line_num = 0
    for line in eachline(open(filename))
        line_num += 1
        if ismatch(ignore_rgx, line)
            continue
        end
        # split and choose columns
        line_array =  map(x->strip(x), split(line,'\t') )

        if seq_id_format == "ucsc"
            push!(seq_id, ens2ucsc( line_array[1]) )
        else
           push!(seq_id, line_array[1] )
        end
        push!(positions, parseint( line_array[2]))
        push!(ref,line_array[4][1] )
        push!(variant,line_array[5][1] )
        if (line_num % 1000000) == 0
            Lumberjack.info("read $line_num lines")
        end
    end
    return DataFrame(seq_id = seq_id, position = positions, ref=ref, variant=variant)
end

# substitutes ref allele with variant allele
function substitute_snps!(seq,snps::DataFrame)
    num_substituted=0
    for row=1:nrow(snps)
        pos = snps[:position][row]
        snp_ref_nuc = snps[:ref][row]
        snp_variant_nuc = snps[:variant][row]
        seq_ref_nuc = seq[pos]
        if uppercase(seq_ref_nuc) == uppercase(snp_ref_nuc)
            seq[pos]=snp_variant_nuc
            num_substituted +=1
        else
            Lumberjack.warn("Got inconsistent reference nucleotide between vcf file and ref sequence at position $pos.\nsnp ref nuc:$snp_ref_nuc seq\nseq ref nuc:$seq_ref_nuc")
        end
    end
    return (num_substituted,nrow(snps))
end

function check_variant_alleles(seq,snps::DataFrame)
    num_substituted=0
    num_ref=0
    num_other=0
    nuc_other_type = Dict()
    for row=1:nrow(snps)
        pos = snps[:position][row]
        snp_ref_nuc = snps[:ref][row]
        snp_variant_nuc = snps[:variant][row]
        seq_ref_nuc = seq[pos]
        if uppercase(seq_ref_nuc) == uppercase(snp_ref_nuc)
            num_ref +=1
        elseif uppercase(seq_ref_nuc) == uppercase( snp_variant_nuc )
            num_substituted +=1
        else
            num_other +=1
            if haskey(nuc_other_type, seq_ref_nuc)
                nuc_other_type[seq_ref_nuc] +=1
            else
                nuc_other_type[seq_ref_nuc]=1
            end
       end
    end
    return (num_substituted,num_ref,num_other,nuc_other_type)
end

function check_variant_alleles_in_fasta_file(fastafile,vcf_filename; seq_id_type="ucsc" )
    Lumberjack.info("Start reading VCF file")
    snp_df = read_snp_positions_from_columns(vcf_filename,seq_id_type)
    fr = FastaReader{Vector{Char}}(fastafile)
    Lumberjack.info("Start checking snp substitution of fasta file")
    for (desc, seq) in fr
         seqlength = length(seq)
         Lumberjack.info("processing $desc - length: $seqlength")
         seqid_df = snp_df[ snp_df[:seq_id] .== desc,:]
         if nrow(seqid_df) > 0
             (num_substituted,num_ref,num_other,nuc_other_type) = check_variant_alleles(seq, seqid_df )
              nuc_type_nums_str = join( map(x->join(x,":"),nuc_other_type), " ")
              Lumberjack.info("num substituted: $num_substituted , num_ref: $num_ref, num_other: $num_other $nuc_type_nums_str")
         end
    end
    Lumberjack.info("read $(fr.num_parsed) entries")
end

function substitute_snps_in_fasta_file(filename,vcf_filename , outfilename)
    Lumberjack.info("Start reading VCF file")
    snp_df = read_snp_positions_from_columns(vcf_filename)
    fr = FastaReader{Vector{Char}}(filename)
    fw = FastaWriter(outfilename)
    Lumberjack.info("Start snp substitution of fasta file")
    for (desc, seq) in fr
         seqlength = length(seq)
         Lumberjack.info("processing $desc - length: $seqlength")
         seqid_df = snp_df[ snp_df[:seq_id] .== desc,:]
         if nrow(seqid_df) > 0
             (num_substituted,total_snps) = substitute_snps!(seq, seqid_df )
              Lumberjack.info("num substituted: $num_substituted of $total_snps")
         end
         writeentry(fw, desc, seq)
    end
   Lumberjack.info("read $(fr.num_parsed) entries")
   close(fw)
end

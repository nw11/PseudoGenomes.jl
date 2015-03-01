using FastaIO
using DataFrames
using Lumberjack
# Gives back a dictionary of dataframes
function read_vcf_columns(filename::String, cols, sequence_id_col=1 )
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

# See vcf-file1.vcf for the kind of file you would expect
# this function to take.
# This is annoying. we don't have a vcf parser, and
# also dataframes does not allow choosing of columns on readtable
# so we don't get its type inference that would pretty much do the job
function read_snp_positions_from_columns(filename::String, sequence_id_col=1, seq_id_prefix="chr" )
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
        push!(seq_id, seq_id_prefix*line_array[1])
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
       end
    end
    return (num_substituted,num_ref,num_other)
end

function check_variant_alleles_in_fasta_file(filename,vcf_filename )
    Lumberjack.info("Start reading VCF file")
    snp_df = read_snp_positions_from_columns(vcf_filename)
    fr = FastaReader{Vector{Char}}(filename)
    Lumberjack.info("Start snp substitution of fasta file")
    for (desc, seq) in fr
         seqlength = length(seq)
         Lumberjack.info("processing $desc - length: $seqlength")
         seqid_df = snp_df[ snp_df[:seq_id] .== desc,:]
         if nrow(seqid_df) > 0
             (num_substituted,num_ref,num_other) = check_variant_alleles(seq, seqid_df )
              Lumberjack.info("num substituted: $num_substituted , num_ref: $num_ref, num_other: $num_other")
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

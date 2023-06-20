
# Counting Codons for filepaths
const TABLE = let
    table = fill(0xff, 2^8)
    for (n, v) in [('A', 0), ('C', 1), ('G', 2), ('T', 3), ('U', 3)]
        table[UInt8(n)+1] = v
        table[UInt8(lowercase(n))+1] = v
    end
    for n in "SWKYMRBDHVN"
        table[UInt8(n)+1] = 0xf0
        table[UInt8(lowercase(n))+1] = 0xf0
    end
    Tuple(table)
end

function count_codons!(vector::AbstractVector{<:Integer}, seq::AbstractString, rem::Integer)
    fill!(vector, 0)
    mask = UInt(1 << 6 - 1)
    remaining = rem
    kmer = UInt8(0)
    for codeunit in codeunits(seq)
        value = TABLE[codeunit+0x01]
        if value == 0xff
            throw(DomainError(codeunit, "Cannot interpret as nucleotide"))
        elseif value == 0xf0
            remaining += 3 # This will "skip" this Kmer entirely, but keep us in frame
        else
            remaining -= 1
            kmer = (kmer << 2 | value) & mask # This needed to be moved outside of the second if block
        end
        if remaining < 1
            @inbounds vector[kmer+1] += 1
            remaining = 3 # This needed to be reset per codon "cycle"
        end
    end
end

"""
    count_codons(path::AbstractString, remove_start::Bool = false, threshold::Integer = 0)
    count_codons(stream::IO, remove_start::Bool = false, threshold::Integer = 0)
    count_codons(reader::FASTAReader, remove_start::Bool = false, threshold::Integer = 0)
    count_codons(sequences::Vector{<:NucSeq}, names::Vector{String} = String[], remove_start::Bool = false, threshold::Integer = 0)
    count_codons(sequence::NucSeq, remove_start::Bool = false, threshold::Integer = 0)
Read a fasta file or BioSequence and return the occurence of each codon for each gene or sequence.

# Arguments
- `path` or `stream` or `reader` or `sequence(s)`: Fasta sequence to analyze. This can be a path to a fasta file of sequences, an IOStream, an open FASTAReader, or a BioSequences nucleotide sequence, or a vector of nucleotide sequences. Note that count_codons isn't identifying ORFs - make sure these are actual CDSs in frame.
- `remove_start`: Whether to ignore the initial start codon
- `threshold`: Minimum length of the sequence *in codons* to be returned in the results.

# Output
If providing a single sequence, the result will be a 64x1 Matrix, which corresponds to the 64 codons in alphabetical order. If you want a list of the codons in alphabetical order, this is stored in `CUBScout.DEFAULT_CodonDict.codons`. If analyzing a fasta file or a vector of sequences, the result will be a tuple. The first element of the tuple is a 64xn matrix, where n = # of sequences above the threshold. The second element is a list of corresponding names for each column. The third element is a Boolean vector where `true` corresponds to sequences which did pass the threshold, and `false` is sequences which did not pass the threshold and so are not included in the results matrix. Names are pulled from fasta files and IO streams by default; if you would like to provide a vector of IDs or names when providing a `Vector{<:NucSeq}`, you can.

# Examples
```jldoctest
julia> example_dna = dna"ATGAAAATGAACTTTTGA"

julia> count_codons(example_dna) |> first
1

julia> result = count_codons(EXAMPLE_DATA_PATH);

julia> first(result[1], 5)
5-element Vector{Int32}:
 32
  7
  6
 14
 11
```
"""
function count_codons(reader::FASTAReader, remove_start::Bool = false, threshold::Integer = 0)
    buffer = zeros(Int, 64)
    result = Int32[]
    names = String[]
    length_passes = Bool[]
    rem = remove_start ? 6 : 3
    for record in reader
        count_codons!(buffer, sequence(record), rem)
        length_pass = sum(buffer) > threshold
        push!(length_passes, length_pass)
        if length_pass
            @inbounds append!(result, buffer)
            @inbounds push!(names, identifier(record))
        end
    end
    @inbounds (reshape(result, 64, :), names, length_passes)
end

function count_codons(path::AbstractString, remove_start::Bool = false, threshold::Integer = 0)
    open(FASTAReader, path; copy = false) do reader
        count_codons(reader, remove_start, threshold)
    end
end

function count_codons(stream::IO, remove_start::Bool = false, threshold::Integer = 0)
    FASTAReader(stream) do reader
        count_codons(reader, remove_start, threshold)
    end
end


# Counting codons for BioSequences
function count_codons(sequence::NucSeq, remove_start::Bool = false)
    cod_space = zeros(Int, (4,4,4)) 
    remaining = remove_start ? 6 : 3
    index = Int8[0,0,0]
    for nuc in sequence
        if remaining > 3 
            remaining -= 1
            continue
        end
        bit = BioSequences.encoded_data(nuc)
        if count_ones(bit) == 1 # 34 here
              @inbounds index[remaining] = trailing_zeros(bit) + 1 
              remaining -= 1
        else
         remaining -= 3
        end
        if remaining < 1 #Another 23
            @inbounds @views cod_space[index[1], index[2], index[3]] += 1 
            remaining = 3
        end
    end
    return @inbounds @views reshape(cod_space, 64, 1)
end

function count_codons!(cod_array::AbstractArray{<:Integer}, index::AbstractArray{<:Integer}, seq::NucSeq, rem::Integer)
    fill!(cod_array, 0)
    remaining = rem
    for nuc in seq
        if remaining > 3 
            remaining -= 1
            continue
        end
        bit = BioSequences.encoded_data(nuc)
        if count_ones(bit) == 1 # 34 here
              @inbounds index[remaining] = trailing_zeros(bit) + 1 
              remaining -= 1
        else
         remaining -= 3
        end
        if remaining < 1 #Another 23
            @inbounds @views cod_array[index[1], index[2], index[3]] += 1 
            remaining = 3
        end
    end
    return @inbounds @views reshape(cod_array, 64, 1)
end



function count_codons(sequences::Vector{<:NucSeq}; names::Union{Vector{String}, Nothing} = nothing, remove_start::Bool = false, threshold::Integer = 0)
    buffer = zeros(Int, (4,4,4))
    i_array = zeros(Int, 3)
    result = Int32[]
    length_passes = Bool[]
    rem = remove_start ? 6 : 3
    for cds in sequences
        count_codons!(buffer, i_array, cds, rem)
        length_pass = sum(buffer) > threshold
        push!(length_passes, length_pass)
        if length_pass
            @inbounds append!(result, buffer)
        end
    end
    name_vec = @inbounds @views isnothing(names) ? nothing : names[length_passes]
    @inbounds (reshape(result, 64, :), name_vec, length_passes)
end

function countsbyAA(count_matrix, dict_uniqueI)
    aa_matrix = Matrix{Int64}(undef, length(dict_uniqueI), size(count_matrix, 2))
    for (i, aa) in enumerate(dict_uniqueI)
        for (j, col) in enumerate(eachcol(selectdim(count_matrix, 1, aa)))
            @inbounds aa_matrix[i, j] = sum(col)
        end
    end
    return aa_matrix
end


function normFrequency(count_matrix, AAcount_matrix, seq_length::Integer, dict_uniqueI)
    freq_matrix = zeros(Float64, size(count_matrix, 1), seq_length)
    AAs = length(dict_uniqueI)
    for column = 1:seq_length
        for row = 1:AAs
            @inbounds freq_matrix[dict_uniqueI[row], column] =
                @views count_matrix[dict_uniqueI[row], column] / AAcount_matrix[row, column]
        end
    end
    return freq_matrix
end

function normTotalFreq(count_matrix, AAcount_matrix, dict_uniqueI)
    rowsums_codon = @views sum(count_matrix, dims = 2)
    rowsums_aa = @views sum(AAcount_matrix, dims = 2)
    freq_vector = zeros(Float64, size(count_matrix, 1))
    for (j, aacount) in enumerate(rowsums_aa)
        @inbounds freq_vector[dict_uniqueI[j]] =
            @views rowsums_codon[dict_uniqueI[j], :] / (aacount)
    end
    return freq_vector
end

function scuo_freq(count_matrix, AA_count_matrix, seq_length, dict_uniqueI)
    freq_matrix = zeros(Float64, size(AA_count_matrix, 1), seq_length)
    AAs = length(dict_uniqueI)
    for column = 1:seq_length
        for row = 1:AAs
            @inbounds freqs =
                count_matrix[dict_uniqueI[row], column] / AA_count_matrix[row, column]
            @inbounds vals = @. -freqs * log10(freqs)
            vals = map((x) -> isnan(x) ? 0.0 : x, vals)
            @inbounds freq_matrix[row, column] = sum(vals)

        end
    end
    return freq_matrix
end

function correction_term(AAcount_matrix, length_vector, dict_deg)
    cor = Float64[]
    for (seq, len) in zip(eachcol(AAcount_matrix), length_vector)
        @inbounds push!(cor, (sum((seq .> 0) .* (dict_deg .- 1)) / len) - 0.5)
    end
    return cor
end

function enc_pi(count_matrix, AAcount_matrix, seq_length::Integer, dict_uniqueI)
    enc_pi = zeros(Float64, length(dict_uniqueI), seq_length)
    AAs = length(dict_uniqueI)
    for column = 1:seq_length
        for row = 1:AAs
            @inbounds freq =
                count_matrix[dict_uniqueI[row], column] / AAcount_matrix[row, column]
            @inbounds enc_pi[row, column] = sum(freq .^ 2)
        end
    end
    return enc_pi
end


# Working on effNC
function eFFNc(fa_matrix, dict_deg)
    @inbounds red = unique(dict_deg)[unique(dict_deg).!=1]
    avgs = map(red) do x
        rows = findall(y -> x .== y, dict_deg)
        @inbounds avg =
            sum(fa_matrix[rows, :], dims = 1) ./ sum(fa_matrix[rows, :] .!= 0, dims = 1)
        @inbounds x == 3 || (avg[avg.==0] .= (1 / x))
        return length(rows) ./ avg
    end
    avgs = reduce(vcat, avgs)
    if any(avgs .== 0)
        threes = findfirst(x -> x == 3, red)
        twos = findfirst(x -> x == 2, red)
        fours = findfirst(x -> x == 4, red)
        cols = findall(x -> x .== 0, avgs[threes, :])
        @inbounds avgs[threes, cols] .=
            (
                avgs[twos, cols] / sum(dict_deg .== 2) +
                avgs[fours, cols] / sum(dict_deg .== 4)
            ) / 2
    end
    enc = sum(dict_deg .== 1) .+ sum(avgs, dims = 1)
    @inbounds enc[enc.>61] .= 61

    return enc
end

"""
    find_seqs(path::AbstractString, match_pattern::Regex)
Read a fasta file at `path` and query the *description* field for a given Regex `match_pattern`. These results can be supplied in either the reference tuples (for codon usage bias functions) or reference vectors (for expressivity measures).

# Examples
```jldoctest
julia> find_seqs(EXAMPLE_DATA_PATH, r"ribosomal")[1:5]
5-element Vector{Bool}:
 0
 0
 0
 0
 0
 ```
"""
function find_seqs(path::AbstractString, match_pattern::Regex)
    open(FASTAReader, path; copy = false) do reader
        match_vector = Bool[]
        for record in reader
            @inbounds push!(match_vector, occursin(match_pattern, description(record)))
        end
        return match_vector
    end
end

"""
    seq_names(path::AbstractString)
Read a fasta file at `path` and return the *name* fields. Just adds convenience on top of FASTX functions.

# Examples
```jldoctest
julia> seq_name_vector = seq_names(EXAMPLE_DATA_PATH);

julia> seq_name_vector[1]
"lcl|NC_000964.3_cds_NP_387882.1_1"
```
"""
function seq_names(path::AbstractString)
    open(FASTAReader, path; copy = false) do reader
        name_vector = String[]
        for record in reader
            @inbounds push!(name_vector, identifier(record))
        end
        return name_vector
    end
end

"""
    seq_descriptions(path::AbstractString)
Read a fasta file at `path` and return the *description* fields. Just adds convenience on top of FASTX functions.

# Examples
```jldoctest
julia> seq_descr = seq_descriptions(EXAMPLE_DATA_PATH);

julia> seq_descr[1]
"lcl|NC_000964.3_cds_NP_387882.1_1 [gene=dnaA] [locus_tag=BSU_00010] [db_xref=EnsemblGenomes-Gn:BSU00010,EnsemblGenomes-Tr:CAB11777,GOA:P05648,InterPro:IPR001957,InterPro:IPR003593,InterPro:IPR010921,InterPro:IPR013159,InterPro:IPR013317,InterPro:IPR018312,InterPro:IPR020591,InterPro:IPR024633,InterPro:IPR027417,PDB:4TPS,SubtiList:BG10065,UniProtKB/Swiss-Prot:P05648] [protein=chromosomal replication initiator informational ATPase] [protein_id=NP_387882.1] [location=410..1750] [gbkey=CDS]"
```
"""
function seq_descriptions(path::AbstractString)
    open(FASTAReader, path; copy = false) do reader
        desc_vector = String[]
        for record in reader
            @inbounds push!(desc_vector, description(record))
        end
        return desc_vector
    end
end

# Functions for dealing with NaNs when necessary
function remove_nan(x, replacement)
    isnan(x) ? replacement : x
end

nanmean(x) = mean(filter(!isnan, x))
nanmean(x, y) = mapslices(nanmean, x, dims = y)

"""
    codon_frequency(codon_counts::Matrix{<:Integer}, form::String, dict::CodonDict = DEFAULT_CodonDict)
Calculate codon frequency from a matrix of codon counts. Accepts as its first argument a `Matrix{<:Integer}` which is a product of `count_codons()`. `form` can be one of four options:
-`net_genomic`: Frequency of each codon across entire genome (matrix).
-`net_gene`: Frequency of each codon within each gene (column).
-`byAA_genomic`: Frequency of each codon within each amino acid across the entire genome (matrix).
-`byAA_gene`: Frequency of each codon within each amino acid within each gene (column). 

If using an alternative genetic code, a custom `CodonDict` can be provided.

# Examples
```jldoctest
julia> codon_counts = count_codons(EXAMPLE_DATA_PATH);

julia> count_matrix = codon_counts[1];

julia> codon_frequency(count_matrix, "net_genomic")[1:5]
5-element Vector{Float64}:
 0.04941242971710299
 0.017114892645228374
 0.021009352696846777
 0.022269444158755328
 0.022257296747490142

julia> codon_frequency(count_matrix, "net_genomic") |> size
(64, 1)

julia> codon_frequency(count_matrix, "net_gene") |> size
(64, 4237)
```
"""
function codon_frequency(codon_counts::Matrix{<:Integer}, form::String, dict::CodonDict = DEFAULT_CodonDict)
    form in ("net_genomic", "byAA_genomic", "net_gene", "byAA_gene") || error("""Invalid form. Please provide. Acceptable forms include "net_genomic", "byAA_genomic", "net_gene", or "byAA_gene".""")
    if form == "net_genomic"
        return @views sum(codon_counts, dims = 2) ./ sum(codon_counts)
    elseif form == "net_gene"
        geneNet = Float64[]
        lengths = @views transpose(sum(codon_counts, dims = 1))
        for (x, l) in zip(eachcol(codon_counts), lengths)
            append!(geneNet, x ./ l)
        end
        return reshape(geneNet, 64, :)
    else
        countAA = countsbyAA(codon_counts, dict.uniqueI)
        if form == "byAA_genomic"
            freq = normTotalFreq(codon_counts, countAA, dict.uniqueI)
            return @views remove_nan.(freq, 0)
        end
        seqs = @views size(codon_counts, 2) 
        freq = normFrequency(codon_counts, countAA, seqs, dict.uniqueI)
        return @views remove_nan.(freq, 0)
    end

end
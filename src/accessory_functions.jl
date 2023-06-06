
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

function count_codons!(vector::AbstractVector{<:Integer}, seq::AbstractString, rem)
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


function count_codons(reader::FASTAReader, remove_start, threshold)
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

function count_codons(path::AbstractString, remove_start, threshold)
    open(FASTAReader, path; copy = false) do reader
        count_codons(reader, remove_start, threshold)
    end
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

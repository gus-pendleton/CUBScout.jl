# B
"""
    b(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
    b(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
Calculate B from Karlin and Mrazek, 1996. 

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple sets of sequences and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of sequences.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = b(EXAMPLE_DATA_PATH); # Calculate measure on example dataset

julia> result_300 = b(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> round.(result.self[1:5], digits = 6)
5-element Vector{Float64}:
 0.209127
 0.328976
 0.223653
 0.539114
 0.249196

julia> b(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> b(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons

julia> all_genes = find_seqs(EXAMPLE_DATA_PATH, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> b(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,)); # Calculate using ribosomal genes as a reference subset

julia> b(EXAMPLE_DATA_PATH, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,)); # Calculate using all genes and ribosomal genes as a reference subset
```
"""
function b(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true, 64)
    end

    return b(sequences, ref_seqs, uniqueI, stop_mask, rm_start, threshold, names)
end

function b(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true, 64)
    end
    if isempty(ref_seqs) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                b(sequences[i], ref_seqs, uniqueI, stop_mask, rm_start, threshold, names)
        end
    elseif isempty(ref_seqs)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                b(sequences[i], ref_seqs, uniqueI, stop_mask, rm_start, threshold, names[i])
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                b(sequences[i], ref_seqs[i], uniqueI, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                b(sequences[i], ref_seqs[i], uniqueI, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function b(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_seqs,
    dict_uniqueI::Vector{Vector{Int32}},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end# Count codons in each gene 
    @inbounds count_matrix = @views counts[1] # Count matrix 64 (codons) x n sequences
    @inbounds names = @views counts[2] # Names of each fasta sequence
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    lengths = @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),))
    else
        @inbounds ref_seqs = @views map(x -> x[counts[3]], ref_seqs)
    end
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Sum total codons for each amino acid for each sequence
    @inbounds pa = @views map(x -> x ./ lengths, eachrow(countAA)) # Find frequency of AA in each gene
    pa = transpose(reduce(hcat, pa))
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Calculate frequency of each codon for each amino acid in each sequence
    @inbounds normsetfreqs = @views map(
        x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI),
        ref_seqs,
    ) # Calculate frequency of each codon for each amino acid in reference subset
    @inbounds dts = @views map(x -> abs.((normfreq .- x)), normsetfreqs) # Subtract the reference frequency of codon from the frequency of that codon within each gene
    dts = map(dts) do y
        map((x) -> remove_nan(x, 0.0), y) # Replace nans with 0s (will be summed later)
    end
    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = @views sum(dt[aa, :], dims = 1) # Sum up contribution of dt for each amino acid
            @inbounds ba[i, :] = row
        end
        ba
    end
    bs = map(bas) do ba
        @inbounds vec(sum(ba .* pa, dims = 1)) # Multiply ba by pa and sum for each gene sequence
    end
    return (bs..., Identifier = names)
end

# ENC
"""
    enc(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    enc(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate ENC from Wright, 1990.  

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = enc(EXAMPLE_DATA_PATH); # Run ENC on example dataset

julia> round.(result.ENC[1:5], digits = 6)
5-element Vector{Float64}:
 56.787282
 52.725947
 59.287949
 52.296686
 55.262981

julia> result_300 = enc(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.ENC)
3801

julia> length(result_300.ENC)
1650

julia> enc(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> enc(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons
```
"""
function enc(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return enc(sequences, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function enc(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc(sequences[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc(sequences[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function enc(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1] # This returns the codon count matrix
    @inbounds names = @views counts[2] # This is the names for each sequence in the file
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Sum total codons for each amino acid for each sequence 
    pi_vec = enc_pi(count_matrix, countAA, seqs, dict_uniqueI) # Calculate pi statistic for each gene
    @inbounds fa = @views @. (countAA * pi_vec - 1) / (countAA - 1) # Calculate Fa
    fa[isnan.(fa)] .= 0.0 # Replace NaN with 0.0 (okay because will sum next)
    res = vec(eFFNc(fa, dict_deg)) # Calculate Nc for each gene
    return (ENC = res, Identifier = names)
end

# ENC Prime
"""
    enc_p(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
    enc_p(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
Calculate ENC' from Novembre, 2002. 

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple sets of sequences and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of sequences.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = enc_p(EXAMPLE_DATA_PATH); # Calculate measure on example dataset

julia> result_300 = enc_p(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> round.(result.self[1:5], digits = 6)
5-element Vector{Float64}:
 61.0
 59.369798
 60.749462
 61.0
 61.0

julia> enc_p(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> enc_p(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons

julia> all_genes = find_seqs(EXAMPLE_DATA_PATH, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> enc_p(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,)); # Calculate using ribosomal genes as a reference subset

julia> enc_p(EXAMPLE_DATA_PATH, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,)); # Calculate using all genes and ribosomal genes as a reference subset
```
"""
function enc_p(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return enc_p(sequences, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function enc_p(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isempty(ref_seqs) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc_p(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    elseif isempty(ref_seqs)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc_p(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc_p(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                enc_p(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end


function enc_p(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_seqs,
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1] # Count matrix 64 (codons) x n sequences
    @inbounds names = @views counts[2] # Names of each fasta sequence
    @inbounds ount_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),))
    else
        @views ref_seqs = @views map(x -> x[counts[3]], ref_seqs)
    end
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Sum total codons for each amino acid for each sequence
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find frequency of each codon within aino acid for each gene
    @inbounds normsetfreqs = @views map(
        x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI),
        ref_seqs,
    ) # Find the frequency of each codon across reference subset
    @inbounds dts = @views map(x -> (@. (normfreq - x)^2 / x), normsetfreqs) # Calculate deviation from reference set for each codon
    dts = map(dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end # Remove NaNs
    chisums = map(dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = @views sum(dt[aa, :], dims = 1) # Sum up deviations for each amino acid
            @inbounds chisum[i, :] = row
        end
        chisum
    end
    chisqs = map(x -> x .* countAA, chisums) # Calculate chi-squared values
    fas = map(chisqs) do chisq
        @inbounds fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg) # Factor in degeneracy to calculate Fa
        @inbounds fa[countAA.<5] .= 0.0
        fa
    end
    res = map(x -> vec(eFFNc(x, dict_deg)), fas) # Calculate Nc
    return (res..., Identifier = names)
end

# MCB 
"""
    mcb(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
    mcb(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
Calculate MCB from Urrutia and Hurst, 2001.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple sets of sequences and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of sequences.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = mcb(EXAMPLE_DATA_PATH); # Calculate measure on example dataset

julia> result_300 = mcb(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> round.(result.self[1:5], digits = 6)
5-element Vector{Float64}:
 0.087211
 0.178337
 0.189682
 0.24012
 0.149869

julia> mcb(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> mcb(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons

julia> all_genes = find_seqs(EXAMPLE_DATA_PATH, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> mcb(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,)); # Calculate using ribosomal genes as a reference subset

julia> mcb(EXAMPLE_DATA_PATH, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,)); # Calculate using all genes and ribosomal genes as a reference subset
```
"""
function mcb(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return mcb(sequences, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function mcb(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isempty(ref_seqs) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                mcb(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    elseif isempty(ref_seqs)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                mcb(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                mcb(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                mcb(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function mcb(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_seqs,
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1] # Count matrix 64 (codons) x n sequences
    @inbounds names = @views counts[2] # Names of each fasta sequence
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),))
    else
        @inbounds ref_seqs = @views map(x -> x[counts[3]], ref_seqs)
    end
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Sum total codons for each amino acid for each sequence
    @inbounds A = countAA[dict_deg.>1, :] .> 0
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find frequency of each codon within each amino acid for each gene
    @inbounds normsetfreqs = @views map(
        x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI),
        ref_seqs,
    ) # Find the frequency of each codon across reference subset
    @inbounds dts = map(x -> (@. (normfreq - x)^2 / x), normsetfreqs)
    dts = map(dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end # Remove NaNs
    no_counts = count_matrix .<= 0
    dts = map(dts) do dt
        @inbounds dt[no_counts] .= 0
        dt
    end

    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = @views sum(dt[aa, :], dims = 1)
            @inbounds ba[i, :] = row
        end
        ba
    end

    mcbs = map(bas) do ba
        mat1 = @. ba[dict_deg>1, :] * log10(countAA[dict_deg>1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end
    return (mcbs..., Identifier = names)
end

# MILC
"""
    milc(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
    milc(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
Calculate MILC from Supek and Vlahovicek, 2005.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple sets of sequences and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of sequences.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = milc(EXAMPLE_DATA_PATH); # Calculate measure on example dataset

julia> result_300 = milc(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> round.(result.self[1:5], digits = 6)
5-element Vector{Float64}:
 0.494826
 0.583944
 0.499472
 0.635493
 0.543935

julia> milc(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> milc(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons

julia> all_genes = find_seqs(EXAMPLE_DATA_PATH, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> milc(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,)); # Calculate using ribosomal genes as a reference subset

julia> milc(EXAMPLE_DATA_PATH, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,)); # Calculate using all genes and ribosomal genes as a reference subset
```
"""
function milc(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return milc(sequences, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function milc(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isempty(ref_seqs) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                milc(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    elseif isempty(ref_seqs)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                milc(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                milc(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                milc(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function milc(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_seqs,
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1] # Count matrix 64 (codons) x n sequences
    @inbounds names = @views counts[2] # Names of each fasta sequence
    @inbounds count_matrix = count_matrix[stop_mask, :] # Remove entries if removing stop codons
    lengths = @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),))
    else
        ref_seqs = @views map(x -> x[counts[3]], ref_seqs)
    end
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Sum total codons for each amino acid for each sequence
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find the frequency of each codon for each gene
    @inbounds normsetfreqs =
        map(x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    @inbounds per_codon_mas = map(x -> (@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
        map((x) -> remove_nan(x, 0.0), y)
    end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds @views row = 2 .* sum(pcmas[aa, :], dims = 1) # Calculate ma for each amino acid
            @inbounds ma[i, :] = row
        end
        ma
    end
    milcs = map(mas) do ma
        @inbounds @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    return (milcs..., Identifier = names)
end

# SCUO
"""
    scuo(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    scuo(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate SCUO from Wan et al., 2004.  

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> result = scuo(EXAMPLE_DATA_PATH); # Run SCUO on example dataset

julia> round.(result.SCUO[1:5], digits = 6)
5-element Vector{Float64}:
 0.143121
 0.191237
 0.096324
 0.345211
 0.105744

julia> result_300 = scuo(EXAMPLE_DATA_PATH, threshold = 300); # Increase threshold length

julia> length(result.SCUO)
3801

julia> length(result_300.SCUO)
1650

julia> scuo(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> scuo(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons
```
"""
function scuo(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return scuo(sequences, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function scuo(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                scuo(sequences[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                scuo(sequences[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function scuo(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds @views count_matrix = counts[1] # This is our codon count matrix
    @inbounds @views names = counts[2] # These are the names of each sequence
    @inbounds @views count_matrix = count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Find amino acid count matrix

    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI) # Calculate normalized frequency of each codon 

    Hmax = log10.(dict_deg)

    Oa = @views map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    @inbounds Oa = reduce(hcat, Oa)

    @inbounds Fa = @views countAA ./ sum(countAA[dict_deg.>1, :], dims = 1)

    @inbounds mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    res = vec(sum(mult, dims = 1))
    return (SCUO = res, Identifier = names)
end

"""
    all_cub(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
    all_cub(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80)
Calculate all codon usage bias measures at once. Because many measures require the same initial calculations, this is more efficient than calculating them individually.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> all_cub_results = all_cub(EXAMPLE_DATA_PATH); # Calculate all six codon usage measures on example dataset

julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> all_cub(EXAMPLE_DATA_PATH, ref_seqs = (ribosomal = ribosomal_genes,)); # Calculate all measures using ribosomal genes as a reference subset
```
"""
function all_cub(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    return all_cub(sequences, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
end

function all_cub(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_seqs = (),
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    len = length(sequences)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
    end
    if isempty(ref_seqs) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                all_cub(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    elseif isempty(ref_seqs)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                all_cub(sequences[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                all_cub(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names)
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] =
                all_cub(sequences[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, names[i])
        end
    end
    return results
end

function all_cub(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_seqs,
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
    counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
        count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
    else count_codons(fasta_seq, rm_start, threshold) 
    end # Count codons in each gene 
    @inbounds @views count_matrix = counts[1] # This is our codon count matrix
    @inbounds @views seqnames = counts[2] # These are our sequence names
    @inbounds @views count_matrix = count_matrix[stop_mask, :] # Remove entries if removing stop codons
    lengths = @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),)) # If no reference subset provided, make a "self" subset
    else
        ref_seqs = @views map(x -> x[counts[3]], ref_seqs)
    end
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Count amino acids in each gene
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Calculate codon frequency within each amino acid within each gene
    @inbounds @views normsetfreqs =
        map(x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI), ref_seqs) # Calculate codon frequency within each amino acid across reference subsets
    # Up to this point, all of the measures should have been the same

    # Calculating B
    @inbounds b_pa = map(x -> x ./ lengths, eachrow(countAA))
    b_pa = transpose(reduce(hcat, b_pa))

    @inbounds b_dts = map(x -> abs.((normfreq .- x)), normsetfreqs) # This is good
    b_dts = map(b_dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end
    b_bas = map(b_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds @views row = @views sum(dt[aa, :], dims = 1)
            @inbounds @views ba[i, :] = row
        end
        ba
    end
    B_result = map(b_bas) do ba
        vec(sum(ba .* b_pa, dims = 1))
    end
    # End calculating B (B is the result)

    # Calculate ENC
    pi_vec = enc_pi(count_matrix, countAA, seqs, dict_uniqueI)
    @inbounds enc_fa = @. (countAA * pi_vec - 1) / (countAA - 1)
    @inbounds enc_fa[isnan.(enc_fa)] .= 0.0
    ENC_result = (ENC = vec(eFFNc(enc_fa, dict_deg)),)
    # End calculating ENC (ENC is the result)

    # Calculate ENC'
    @inbounds encp_dts = map(x -> (@. (normfreq - x)^2 / x), normsetfreqs) # Calculate Ma for each codon
    encp_dts = map(encp_dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end # Remove NaNs
    encp_chisums = map(encp_dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = @views sum(dt[aa, :], dims = 1)
            @inbounds chisum[i, :] = row
        end
        chisum
    end
    encp_chisqs = map(x -> x .* countAA, encp_chisums) # This also looks good
    encp_fas = map(encp_chisqs) do chisq
        @inbounds fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg)
        @inbounds fa[countAA.<5] .= 0.0
        fa
    end
    ENCP_result = map(x -> vec(eFFNc(x, dict_deg)), encp_fas)
    #End calculating ENC' (ENCP is the result)

    #Start calculating MCB (can use ENCP dts)
    no_counts = count_matrix .<= 0
    mcb_dts = map(encp_dts) do dt
        @inbounds dt[no_counts] .= 0
        dt
    end

    mcb_bas = map(mcb_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = @views sum(dt[aa, :], dims = 1)
            @inbounds ba[i, :] = row
        end
        ba
    end
    A = countAA[dict_deg.>1, :] .> 0
    MCB_result = map(mcb_bas) do ba
        @inbounds mat1 = @. ba[dict_deg>1, :] * log10(countAA[dict_deg>1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        @inbounds vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end

    # End calculating MCB (MCB is the result)

    # Start calculating MILC
    cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    @inbounds per_codon_mas = map(x -> (@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
        map((x) -> remove_nan(x, 0.0), y)
    end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i, aa) in enumerate(dict_uniqueI)
            @inbounds row = 2 .* sum(pcmas[aa, :], dims = 1) # Calculate ma for each amino acid
            @inbounds ma[i, :] = row
        end
        ma
    end
    MILC_result = map(mas) do ma
        @inbounds @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    # End calculating MILC (MILC is the result)

    # Start calculating SCUO
    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI)

    Hmax = log10.(dict_deg)

    Oa = map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    Oa = reduce(hcat, Oa)

    @inbounds Fa = countAA ./ sum(countAA[dict_deg.>1, :], dims = 1)

    @inbounds mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    SCUO_result = (SCUO = vec(sum(mult, dims = 1)),)

    return (
        B = B_result,
        ENC = ENC_result,
        ENCP = ENCP_result,
        MCB = MCB_result,
        MILC = MILC_result,
        SCUO = SCUO_result,
        Identifier = seqnames,
    )

end

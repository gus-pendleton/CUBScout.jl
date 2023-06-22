
"""
    melp(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, ref_vector::Vector{Bool}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    melp(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, ref_vectors::Vector{Vector{Bool}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate MELP from Supek and Vlahovicek, 2005.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `ref_vector`: reference subset, which is required for `melp`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = melp(EXAMPLE_DATA_PATH, ribosomal_genes); # Calculate MELP on example dataset

julia> round.(result.MELP[1:5], digits = 6)
5-element Vector{Float64}:
 0.929414
 1.007671
 0.922357
 0.951239
 1.029531

julia> melp(EXAMPLE_DATA_PATH, ribosomal_genes, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> melp(EXAMPLE_DATA_PATH, ribosomal_genes, rm_start = true); # Remove start codons
```
"""
function melp(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    milcs = milc(
        sequences,
        dict,
        ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector),
        rm_start = rm_start,
        rm_stop = rm_stop,
        threshold = threshold,
        names = names
    )
    return (MELP = milcs.self ./ milcs.reference, Identifier = milcs.Identifier)
end

function melp(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    ref_vectors::Vector{Vector{Bool}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    ref_tuples = map(x -> (self = fill(true, length(x)), reference = x), ref_vectors)
    milcs = milc(
        sequences,
        dict,
        ref_seqs = ref_tuples,
        rm_start = rm_start,
        rm_stop = rm_stop,
        threshold = threshold,
        names = names
    )
    return map(x -> (MELP = x.self ./ x.reference, x.Identifier), milcs)
end

"""
    e(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, ref_vector::Vector{Bool}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    e(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, ref_vectors::Vector{Vector{Bool}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate E from Karlin and Mrazek, 1996.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `ref_vector`: reference subset, which is required for `e`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = e(EXAMPLE_DATA_PATH, ribosomal_genes); # Calculate E on example dataset

julia> round.(result.E[1:5], digits = 6)
5-element Vector{Float64}:
 0.762317
 1.025839
 0.875954
 0.986498
 1.111275

julia> e(EXAMPLE_DATA_PATH, ribosomal_genes, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> e(EXAMPLE_DATA_PATH, ribosomal_genes, rm_start = true); # Remove start codons
```
"""
function e(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    bs = b(
        sequences,
        dict,
        ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector),
        rm_start = rm_start,
        rm_stop = rm_stop,
        threshold = threshold,
        names = names
    )
    return (E = bs.self ./ bs.reference, Identifier = bs.Identifier)
end

function e(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    ref_vectors::Vector{Vector{Bool}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    rm_start = false,
    rm_stop = false,
    threshold = 80,
)
    ref_tuples = map(x -> (self = fill(true, length(x)), reference = x), ref_vectors)
    bs = b(
        sequences,
        dict,
        ref_seqs = ref_tuples,
        rm_start = rm_start,
        rm_stop = rm_stop,
        threshold = threshold,
        names = names
    )
    return map(x -> (E = x.self ./ x.reference, x.Identifier), bs)
end

"""
    cai(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, ref_vector::Vector{Bool}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    cai(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, ref_vectors::Vector{Vector{Bool}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate CAI from Sharp and Li, 1987.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `ref_vector`: reference subset, which is required for `cai`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = cai(EXAMPLE_DATA_PATH, ribosomal_genes); # Calculate CAI on example dataset

julia> round.(result.CAI[1:5], digits = 6)
5-element Vector{Float64}:
 0.844967
 0.88548
 0.817348
 1.072675
 0.834179

julia> cai(EXAMPLE_DATA_PATH, ribosomal_genes, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> cai(EXAMPLE_DATA_PATH, ribosomal_genes, rm_start = true); # Remove start codons
```
"""
function cai(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
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
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
        aa_names = dict.AA
    end
    return cai(sequences, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, names)
end

function cai(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    ref_vectors::Vector{Vector{Bool}},
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
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
        aa_names = dict.AA
    end

    if isnothing(names)
        Threads.@threads for i = 1:len
        @inbounds results[i] = cai(
            sequences[i],
            ref_vectors[i],
            uniqueI,
            deg,
            stop_mask,
            aa_names,
            rm_start,
            threshold,
            names
        )
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] = cai(
                sequences[i],
                ref_vectors[i],
                uniqueI,
                deg,
                stop_mask,
                aa_names,
                rm_start,
                threshold,
                names[i]
            )
        end
    end
    return results
end


function cai(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    aa_names::Vector{String},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
    count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
else count_codons(fasta_seq, rm_start, threshold) 
end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1]
    @inbounds names = @views counts[2]
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have

    @inbounds ref_seqs = (self = fill(true, seqs), reference = ref_vector[counts[3]])

    countAA = countsbyAA(count_matrix, dict_uniqueI) # This is the same for all measures
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI)
    @inbounds normsetfreqs =  @views map(
        x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI),
        ref_seqs,
    )
    max_aa = fill(0.0, length(aa_names))
    map(dict_uniqueI) do aa
        @inbounds max_aa[aa] .= maximum(normsetfreqs.reference[aa])
    end
    @inbounds max_aa[max_aa.==0] .= 0.5
    @inbounds nodeg = dict_uniqueI[dict_deg.==1]
    @inbounds map(x -> count_matrix[x, :] .= 0, nodeg)
    @inbounds mult = @. log(normfreq / max_aa) * count_matrix
    mult = remove_nan.(mult, 0)
    @inbounds cai_result = vec(exp.(sum(mult, dims = 1) ./ sum(count_matrix, dims = 1)))
    return (CAI = cai_result, Identifier = names)
end


"""
    fop(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, ref_vector::Vector{Bool}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
    fop(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, ref_vectors::Vector{Vector{Bool}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, rm_start = false, rm_stop = false, threshold = 80)
Calculate FOP from Ikemura, 1981.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of sequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `ref_vector`: reference subset, which is required for `fop`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = fop(EXAMPLE_DATA_PATH, ribosomal_genes); # Calculate CAI on example dataset

julia> round.(result.FOP[1:5], digits = 6)
5-element Vector{Float64}:
 0.567816
 0.566845
 0.509695
 0.725
 0.653784

julia> fop(EXAMPLE_DATA_PATH, ribosomal_genes, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> fop(EXAMPLE_DATA_PATH, ribosomal_genes, rm_start = true); # Remove start codons
```
"""
function fop(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
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
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
        aa_names = dict.AA
    end
    return fop(sequences, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, names)
end

function fop(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    ref_vectors::Vector{Vector{Bool}},
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
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true, 64)
        aa_names = dict.AA
    end

    if isnothing(names)
        Threads.@threads for i = 1:len
        @inbounds results[i] = fop(
            sequences[i],
            ref_vectors[i],
            uniqueI,
            deg,
            stop_mask,
            aa_names,
            rm_start,
            threshold,
            names
        )
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] = fop(
                sequences[i],
                ref_vectors[i],
                uniqueI,
                deg,
                stop_mask,
                aa_names,
                rm_start,
                threshold,
                names[i]
            )
        end
    end
    return results
end


function fop(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    ref_vector::Vector{Bool},
    dict_uniqueI::Vector{Vector{Int32}},
    dict_deg::Vector{<:Integer},
    stop_mask::Vector{Bool},
    aa_names::Vector{String},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
    count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
else count_codons(fasta_seq, rm_start, threshold) 
end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1]
    @inbounds names = @views counts[2]
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have

    @inbounds ref_seqs = (self = fill(true, seqs), reference = ref_vector[counts[3]])

    countAA = countsbyAA(count_matrix, dict_uniqueI) # This is the same for all measures
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI)
    @inbounds normsetfreqs = @views map(
        x -> normTotalFreq(count_matrix[:, x], countAA[:, x], dict_uniqueI),
        ref_seqs,
    )
    max_aa = fill(0.0, length(aa_names))
    map(dict_uniqueI) do aa
        @inbounds max_aa[aa] .= maximum(normsetfreqs.reference[aa])
    end
    @inbounds max_aa[max_aa.==0] .= 0.5
    @inbounds nodeg = dict_uniqueI[dict_deg.==1]
    @inbounds map(x -> count_matrix[x, :] .= 0, nodeg)
    @inbounds ra = normfreq ./ max_aa
    count2 = copy(count_matrix)
    @inbounds count2[ra.<0.9] .= 0
    fops = vec(sum(count2, dims = 1) ./ sum(count_matrix, dims = 1))
    return (FOP = fops, Identifier = names)

end

"""
    gcb(sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{String}, Nothing} = nothing, ref_vector = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80)
    gcb(sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}}, dict::CodonDict = DEFAULT_CodonDict; names::Union{Vector{Vector{String}}, Nothing} = nothing, ref_vectors = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80)
Calculate GCB from Merkl, 2003.

# Arguments
- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment. There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions. If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of BioSequences, with each vector corresponding to a genome. `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each fasta file or vector of BioSequences. As such, providing a vector of filepaths (or `Vector{<:Vector{<:NucSeq}}`) as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed.
- `dict`: codon dictionary of type `CodonDict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_CodonDict`
- `names`: An optional vector of names for each sequence. Only relevant if providing a vector of BioSequences, as names are automatically pulled from fasta files. If `sequences` is of type `Vector{<:Vector{<:NucSeq}}`, `names` should be of type `Vector{Vector{String}}`
- `ref_vector`: optional reference subset; by default gcb begins calculations using all genes as a seed. If you want to provide a custom reference set, it should be a vector `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `perc`: percentage of "top hits" which should be used as a reference set in the next iteration. By default set to 0.05. 
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `ALTSTART_CodonDict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = gcb(EXAMPLE_DATA_PATH); # Calculate GCB on example dataset

julia> round.(result.GCB[1:5], digits = 6)
5-element Vector{Float64}:
 -0.058765
 -0.08659
 -0.005496
 -0.065659
 -0.032062

julia> ribo_result = gcb(EXAMPLE_DATA_PATH, ref_vector = ribosomal_genes); # Calculate GCB with ribosomal genes as reference seed example dataset

julia> round.(ribo_result.GCB[1:5], digits = 6)
5-element Vector{Float64}:
 -0.135615
 -0.036687
 -0.169136
 -0.186104
 -0.01653

julia> gcb(EXAMPLE_DATA_PATH, ALTSTART_CodonDict); # Code TTG and CTG as methionine

julia> gcb(EXAMPLE_DATA_PATH, rm_start = true); # Remove start codons
```
"""
function gcb(
    sequences::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{String}, Nothing} = nothing,
    ref_vector = [],
    perc = 0.05,
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
    return gcb(sequences, ref_vector, uniqueI, perc, stop_mask, rm_start, threshold, names)
end

function gcb(
    sequences::Union{Vector{String}, Vector{<:IO}, Vector{<:FASTAReader}, Vector{<:Vector{<:NucSeq}}},
    dict::CodonDict = DEFAULT_CodonDict;
    names::Union{Vector{Vector{String}}, Nothing} = nothing,
    ref_vectors = [],
    perc = 0.05,
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
    if isempty(ref_vectors) & isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] = gcb(
                sequences[i],
                ref_vectors,
                uniqueI,
                perc,
                stop_mask,
                rm_start,
                threshold,
                names,
            )
        end
    elseif isnothing(names)
        Threads.@threads for i = 1:len
            @inbounds results[i] = gcb(
                sequences[i],
                ref_vectors[i],
                uniqueI,
                perc,
                stop_mask,
                rm_start,
                threshold,
                names
            )
        end
    elseif isempty(ref_vectors)
        Threads.@threads for i = 1:len
            @inbounds results[i] = gcb(
                sequences[i],
                ref_vectors,
                uniqueI,
                perc,
                stop_mask,
                rm_start,
                threshold,
                names[i]
            )
        end
    else
        Threads.@threads for i = 1:len
            @inbounds results[i] = gcb(
                sequences[i],
                ref_vectors[i],
                uniqueI,
                perc,
                stop_mask,
                rm_start,
                threshold,
                names[i]
            )
        end
    end
    return results
end


function gcb(
    fasta_seq::Union{String, IO, FASTAReader, Vector{<:NucSeq}},
    refs,
    dict_uniqueI::Vector{Vector{Int32}},
    perc::Real,
    stop_mask::Vector{Bool},
    rm_start::Bool,
    threshold::Integer,
    names::Union{Vector{String}, Nothing},
)
counts = if typeof(fasta_seq) <: Vector{<:NucSeq}
    count_codons(fasta_seq, names = names, remove_start = rm_start, threshold = threshold)
else count_codons(fasta_seq, rm_start, threshold) 
end # Count codons in each gene 
    @inbounds count_matrix = @views counts[1]
    @inbounds names = @views counts[2]
    @inbounds count_matrix = @views count_matrix[stop_mask, :] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    lengths = @views transpose(sum(count_matrix, dims = 1))

    @inbounds seed = isempty(refs) ? fill(true, seqs) : refs[counts[3]] # Make our seed - this will be our initial reference set
    countAA = countsbyAA(count_matrix, dict_uniqueI) # Count 
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI)
    @inbounds normsetfreq =
        @views normTotalFreq(count_matrix[:, seed], countAA[:, seed], dict_uniqueI)

    gcb_prev = fill(0.0, seqs)
    iter = 0
    gcb = []
    diff = false
    # Now we'd enter the repeat loop
    while true
        @inbounds cb = log.(normsetfreq ./ nanmean(normfreq, 2))
        @inbounds cb[normsetfreq.==0] .= -5
        @inbounds gcb = vec(vec(sum(count_matrix .* cb, dims = 1)) ./ lengths)
        diff = all(gcb .== gcb_prev)
        if diff | iter > 6
            break
        end
        iter += 1
        gcb_prev = copy(gcb)
        @inbounds tops = sortperm(gcb, rev = true)[1:convert(Int, trunc(perc * seqs))]
        seed .= false
        @inbounds seed[tops] .= true
        @inbounds normsetfreq =
            @views normTotalFreq(count_matrix[:, seed], countAA[:, seed], dict_uniqueI)
    end
    return (GCB = gcb, Identifier = names)
end

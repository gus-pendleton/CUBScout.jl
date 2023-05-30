# B
"""
    b(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    b(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate B from Karlin and Mrazek, 1996. 

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple filepaths and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of filepaths.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = b(example_data_path) # Run B on example dataset
(self = [0.20912699220973896, 0.3289759448740455, 0.22365336363593893, 0.5391135258658497, 0.24919594143501034, 0.2880358413249049, 0.31200964304415874, 0.34858035204347476, 0.2455189361074733, 0.4690734561271221  …  0.3629137353834403, 0.3621330537227321, 0.4535285720373026, 0.3357858047622507, 0.28183191395624935, 0.2668809561422238, 0.22381338105820905, 0.4034837015709619, 0.3594626865160133, 0.3724863965444541],)

julia> result_300 = b(example_data_path, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> b(example_data_path, altstart_codon_dict) # Use alternative start codons
(self = [0.20897234061622738, 0.33515000264964157, 0.23659038285006437, 0.5444798345895256, 0.2510726777670733, 0.2931440715811394, 0.32097661134289895, 0.35705614480228676, 0.25452296343822073, 0.513313870450466  …  0.34414371567411556, 0.38229515825882665, 0.4592524704597901, 0.3399477982926337, 0.29297757306048133, 0.2680028918895221, 0.2168486105068708, 0.414543030746344, 0.3829702745346273, 0.39870546723886807],)

julia> b(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ self      Identifier                         File                              
      │ Float64   String                             String                            
──────┼────────────────────────────────────────────────────────────────────────────────
    1 │ 0.209127  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.328976  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.223653  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.539114  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮                      ⋮                                  ⋮
 3799 │ 0.403484  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.359463  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.372486  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                      3794 rows omitted

julia> all_genes = find_seqs(example_data_path, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> b(example_data_path, ref_seqs = (ribosomal = ribosomal_genes,), dataframe = true) # Calculate B using ribosomal genes as a reference subset
3801×3 DataFrame
  Row │ ribosomal  Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │  0.274331  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │  0.32069   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │  0.255325  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │  0.546493  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │  0.406673  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │  0.375857  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │  0.437981  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> b(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate B using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self      ribosomal  Identifier                         File                              
      │ Float64   Float64    String                             String                            
──────┼───────────────────────────────────────────────────────────────────────────────────────────
    1 │ 0.209127   0.274331  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.328976   0.32069   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.223653   0.255325  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.539114   0.546493  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮          ⋮                      ⋮                                  ⋮
 3799 │ 0.403484   0.406673  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.359463   0.375857  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.372486   0.437981  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                 3794 rows omitted
```
"""
function b(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end

    return b(filepath, ref_seqs, uniqueI, stop_mask, rm_start, threshold, dataframe)
end

function b(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = b(filepaths[i], ref_seqs, uniqueI, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = b(filepaths[i], ref_seqs[i], uniqueI, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function b(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # Count matrix 64 (codons) x n sequences
    names = counts[2] # Names of each fasta sequence
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),)) 
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Sum total codons for each amino acid for each sequence
    pa = map(x-> x ./ lengths, eachrow(countAA)) # Find frequency of AA in each gene
    pa = transpose(reduce(hcat, pa))
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Calculate frequency of each codon for each amino acid in each sequence
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Calculate frequency of each codon for each amino acid in reference subset
    dts = map(x->abs.((normfreq .- x)), normsetfreqs) # Subtract the reference frequency of codon from the frequency of that codon within each gene
    dts = map(dts) do y
        map((x) -> remove_nan(x, 0.0), y) # Replace nans with 0s (will be summed later)
    end
    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1) # Sum up contribution of dt for each amino acid
            ba[i,:] = row
        end
        ba
    end
    bs = map(bas) do ba
         vec(sum(ba .* pa, dims = 1)) # Multiply ba by pa and sum for each gene sequence
    end
    dataframe && return (df = DataFrame(bs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df) # Format in dataframe output, if dataframe = true
    return bs
end

# ENC
"""
    enc(filepath::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    enc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate ENC from Wright, 1990. 

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = enc(example_data_path) # Run ENC on example dataset
3801-element Vector{Float64}:
 56.787282202547104
 52.725946690067296
 59.287948966886226
 52.29668642771212
 55.26298060679466
[...]

julia> result_300 = enc(example_data_path, threshold = 300); # Increase threshold length

julia> length(result)
3801

julia> length(result_300)
1650

julia> enc(example_data_path, altstart_codon_dict) # Use alternative start codons
3801-element Vector{Float64}:
 57.001523875733255
 52.86344589270121
 59.5010825383314
 52.24194738136219
 55.24666826718284
 [...]

julia> enc(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ ENC      Identifier                         File                              
      │ Float64  String                             String                            
──────┼───────────────────────────────────────────────────────────────────────────────
    1 │ 56.7873  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 52.7259  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 59.2879  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 52.2967  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮                     ⋮                                  ⋮
 3799 │ 55.2294  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 52.584   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 60.1928  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                     3794 rows omitted
```
"""
function enc(filepath::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return enc(filepath, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function enc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    Threads.@threads for i in 1:len
       @inbounds results[i] = enc(filepaths[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    dataframe && return reduce(vcat, results)
    return results
end

function enc(fasta_seq::String, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # This returns the codon count matrix
    names = counts[2] # This is the names for each sequence in the file
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Sum total codons for each amino acid for each sequence 
    pi_vec =  enc_pi(count_matrix, countAA, seqs, dict_uniqueI) # Calculate pi statistic for each gene
    fa =  @. (countAA * pi_vec - 1) / (countAA - 1) # Calculate Fa
    fa[isnan.(fa)] .= 0.0 # Replace NaN with 0.0 (okay because will sum next)
    res = vec(eFFNc(fa, dict_deg)) # Calculate Nc for each gene
    dataframe && return DataFrame(ENC = res, Identifier = names, File = fasta_seq)
    return res
end

# ENC Prime
"""
    enc_p(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    enc_p(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate ENC' from Novembre, 2002.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple filepaths and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of filepaths.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = enc_p(example_data_path) # Run ENC' on example dataset
(self = [61.0, 59.36979815371983, 60.7494622549966, 61.0, 61.0, 56.353402323266224, 55.025304341802055, 57.30607996896261, 61.0, 49.80272180663614  …  58.275539315744126, 61.0, 58.28913535421893, 61.0, 61.0, 59.45750784609185, 61.0, 59.947884174402645, 59.43051836466144, 61.0],)

julia> result_300 = enc_p(example_data_path, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> enc_p(example_data_path, altstart_codon_dict) # Use alternative start codons
(self = [61.0, 59.25393167941246, 60.451328058538294, 61.0, 61.0, 56.577935943942954, 55.48412223261298, 57.35814028261265, 61.0, 48.78039661159036  …  58.93889392717741, 61.0, 59.1935848198233, 61.0, 61.0, 59.587859071790746, 61.0, 58.81213890142517, 59.40548501932952, 61.0],)

julia> enc_p(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ self     Identifier                         File                              
      │ Float64  String                             String                            
──────┼───────────────────────────────────────────────────────────────────────────────
    1 │ 61.0     lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 59.3698  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 60.7495  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 61.0     lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮                     ⋮                                  ⋮
 3799 │ 59.9479  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 59.4305  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 61.0     lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                     3794 rows omitted

julia> all_genes = find_seqs(example_data_path, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> enc_p(example_data_path, ref_seqs = (ribosomal = ribosomal_genes,), dataframe = true) # Calculate ENC' using ribosomal genes as a reference subset
3801×3 DataFrame
  Row │ ribosomal  Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │   61.0     lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │   58.8882  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │   56.4104  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │   61.0     lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │   56.5325  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │   55.6687  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │   61.0     lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> enc_p(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate ENC' using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self     ribosomal  Identifier                         File                              
      │ Float64  Float64    String                             String                            
──────┼──────────────────────────────────────────────────────────────────────────────────────────
    1 │ 61.0       61.0     lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 59.3698    58.8882  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 60.7495    56.4104  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 61.0       61.0     lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮         ⋮                      ⋮                                  ⋮
 3799 │ 59.9479    56.5325  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 59.4305    55.6687  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 61.0       61.0     lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                3794 rows omitted
                                                                       
```
"""
function enc_p(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return enc_p(filepath, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function enc_p(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = enc_p(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = enc_p(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end


function enc_p(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # Count matrix 64 (codons) x n sequences
    names = counts[2] # Names of each fasta sequence
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),)) 
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Sum total codons for each amino acid for each sequence
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find frequency of each codon within aino acid for each gene
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs) # Calculate deviation from reference set for each codon
    dts = map(dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                    end # Remove NaNs
    chisums = map(dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1) # Sum up deviations for each amino acid
            chisum[i,:] = row
        end
        chisum
    end
    chisqs = map(x-> x .* countAA, chisums) # Calculate chi-squared values
    fas = map(chisqs) do chisq
            fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg) # Factor in degeneracy to calculate Fa
            fa[countAA .< 5] .= 0.0
            fa
        end
        res = map(x->vec(eFFNc(x, dict_deg)), fas) # Calculate Nc
    dataframe && return (df = DataFrame(res); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df) # Return dataframe is dataframe = true
    return res
end

# MCB 
"""
    mcb(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    mcb(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate MCB from Urutia and Hurst, 2001.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple filepaths and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of filepaths.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = mcb(example_data_path) # Run MCB on example dataset
(self = [0.08721123763355611, 0.17833660522736233, 0.18968222510423902, 0.2401197552240717, 0.14986851474968427, 0.2726565434395838, 0.2463076527192269, 0.21311734597718654, 0.12911296594333618, 0.33880929702119683  …  0.20493351356705938, 0.1315657625271143, 0.218615505877176, 0.1589526915862505, 0.11884772940321076, 0.16599279317487328, 0.09899721758865292, 0.20119392008973744, 0.2227777749890099, 0.14183390310488245],)

julia> result_300 = mcb(example_data_path, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> mcb(example_data_path, altstart_codon_dict) # Use alternative start codons
(self = [0.08395614019755322, 0.17886564496310903, 0.1843882189861458, 0.2327059977332034, 0.14375956626276776, 0.26153451676836, 0.2277560915069008, 0.2084721025526067, 0.12748812450452368, 0.35284946520551547  …  0.19412496858111208, 0.15378185701040462, 0.22664543233169301, 0.1572460068164919, 0.12228357240597998, 0.16482459147152181, 0.09433344489896123, 0.2180785328818952, 0.22046000262742788, 0.14934890119286454],)

julia> mcb(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ self       Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │ 0.0872112  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.178337   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.189682   lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.24012    lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │ 0.201194   lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.222778   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.141834   lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> all_genes = find_seqs(example_data_path, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> mcb(example_data_path, ref_seqs = (ribosomal = ribosomal_genes,), dataframe = true) # Calculate MCB using ribosomal genes as a reference subset
3801×3 DataFrame
  Row │ ribosomal  Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │ 0.142922   lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.179081   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.306192   lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.296251   lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │ 0.239881   lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.260092   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.197456   lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> mcb(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate MCB using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self       ribosomal  Identifier                         File                              
      │ Float64    Float64    String                             String                            
──────┼────────────────────────────────────────────────────────────────────────────────────────────
    1 │ 0.0872112  0.142922   lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.178337   0.179081   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.189682   0.306192   lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.24012    0.296251   lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮          ⋮                      ⋮                                  ⋮
 3799 │ 0.201194   0.239881   lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.222778   0.260092   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.141834   0.197456   lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                  3794 rows omitted

julia> mcb(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate MCB using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self       ribosomal  Identifier                         File                              
      │ Float64    Float64    String                             String                            
──────┼────────────────────────────────────────────────────────────────────────────────────────────
    1 │ 0.0872112  0.142922   lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.178337   0.179081   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.189682   0.306192   lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.24012    0.296251   lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮          ⋮                      ⋮                                  ⋮
 3799 │ 0.201194   0.239881   lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.222778   0.260092   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.141834   0.197456   lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                  3794 rows omitted
```
"""
function mcb(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return mcb(filepath, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function mcb(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = mcb(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = mcb(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function mcb(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # Count matrix 64 (codons) x n sequences
    names = counts[2] # Names of each fasta sequence
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),)) 
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Sum total codons for each amino acid for each sequence
    A = countAA[dict_deg .> 1,:] .> 0
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find frequency of each codon within each amino acid for each gene
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs)
    dts = map(dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    no_counts = count_matrix .<= 0 
    dts = map(dts) do dt
        dt[no_counts] .= 0
        dt
    end

    bas = map(dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end

    mcbs = map(bas) do ba
        mat1 = @. ba[dict_deg > 1,:] * log10(countAA[dict_deg > 1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end
    dataframe && return (df = DataFrame(mcbs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return mcbs
end

# MILC
"""
    milc(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    milc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate MILC from Super and Vlahovicek, 2005.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple filepaths and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of filepaths.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = milc(example_data_path) # Run MCB on example dataset
(self = [0.49482573202153163, 0.5839439121281993, 0.49947166558087047, 0.6354929447434434, 0.5439352548027006, 0.6104721251245075, 0.6256398806438782, 0.6228376952086359, 0.5355298113407091, 0.7832276821181443  …  0.5968814155010973, 0.5964500002803941, 0.5930680822246766, 0.5412999510428169, 0.49866919389111675, 0.5830959504630727, 0.5139438478694085, 0.6164434557282711, 0.6018041071661588, 0.48775477465069617],)

julia> result_300 = milc(example_data_path, threshold = 300); # Increase threshold length

julia> length(result.self)
3801

julia> length(result_300.self)
1650

julia> milc(example_data_path, altstart_codon_dict) # Use alternative start codons
(self = [0.49430677065204015, 0.5885230358759501, 0.5025039496927455, 0.6371410353128949, 0.5426521549091665, 0.6108102890984481, 0.6217053122088902, 0.6248550493920512, 0.5375433930570145, 0.822162494240932  …  0.5865793308679926, 0.6188550398687176, 0.5942360103321869, 0.5420673855853781, 0.5054979526357374, 0.5824007342085071, 0.5121161561615022, 0.6308328736610525, 0.6146790189624615, 0.5061508729549385],)

julia> milc(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ self      Identifier                         File                              
      │ Float64   String                             String                            
──────┼────────────────────────────────────────────────────────────────────────────────
    1 │ 0.494826  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.583944  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.499472  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.635493  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮                      ⋮                                  ⋮
 3799 │ 0.616443  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.601804  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.487755  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                      3794 rows omitted


julia> all_genes = find_seqs(example_data_path, r""); # Get a vector which is true for all genes

julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> milc(example_data_path, ref_seqs = (ribosomal = ribosomal_genes,), dataframe = true) # Calculate MILC using ribosomal genes as a reference subset
3801×3 DataFrame
  Row │ ribosomal  Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │  0.532406  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │  0.579498  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │  0.541516  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │  0.668068  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │  0.636678  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │  0.623519  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │  0.551733  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

julia> milc(example_data_path, ref_seqs = (self = all_genes, ribosomal = ribosomal_genes,), dataframe = true) # Calculate MILC using all genes and ribosomal genes as a reference subset
3801×4 DataFrame
  Row │ self      ribosomal  Identifier                         File                              
      │ Float64   Float64    String                             String                            
──────┼───────────────────────────────────────────────────────────────────────────────────────────
    1 │ 0.494826   0.532406  lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.583944   0.579498  lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.499472   0.541516  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.635493   0.668068  lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │    ⋮          ⋮                      ⋮                                  ⋮
 3799 │ 0.616443   0.636678  lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.601804   0.623519  lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.487755   0.551733  lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                                 3794 rows omitted

```
"""
function milc(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return milc(filepath, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function milc(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
    Threads.@threads for i in 1:len
       @inbounds results[i] = milc(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    else
        Threads.@threads for i in 1:len
            @inbounds results[i] = milc(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
         end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function milc(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # Count matrix 64 (codons) x n sequences
    names = counts[2] # Names of each fasta sequence
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs) # If no ref_seqs provided, create a "self" tuple
        (ref_seqs = (self = fill(true, seqs),)) 
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Sum total codons for each amino acid for each sequence
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Find the frequency of each codon for each gene
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Find the frequency of each codon across reference subset
    @inbounds cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    @inbounds per_codon_mas = map(x->(@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i,aa) in enumerate(dict_uniqueI)
            row = 2 .* sum(pcmas[aa,:],dims = 1) # Calculate ma for each amino acid
            ma[i,:] = row
        end
        ma
    end
    milcs = map(mas) do ma
        @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    dataframe && return (df = DataFrame(milcs); df[!,:Identifier] = names; df[!,:File] .= fasta_seq; df)
    return milcs
end

# SCUO
"""
    scuo(filepath::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    scuo(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate SCUO from Wan et al., 2004. 

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> result = scuo(example_data_path) # Run SCUO on example dataset
3801-element Vector{Float64}:
 0.14312092935182216
 0.19123738759808445
 0.09632387849329298
 0.34521070127179804
 0.1057438762588421
 0.12037900551695478
[...]

julia> result_300 = scuo(example_data_path, threshold = 300); # Increase threshold length

julia> length(result)
3801

julia> length(result_300)
1650

julia> scuo(example_data_path, altstart_codon_dict) # Use alternative start codons
3801-element Vector{Float64}:
 0.14286111587263958
 0.19315278493814017
 0.0966128845976179
 0.3473543659821751
 0.10792236840320082
 0.12039525638448735
 [...]

julia> scuo(example_data_path, dataframe = true) # Get output in dataframe format
3801×3 DataFrame
  Row │ SCUO       Identifier                         File                              
      │ Float64    String                             String                            
──────┼─────────────────────────────────────────────────────────────────────────────────
    1 │ 0.143121   lcl|NC_000964.3_cds_NP_387882.1_1  /Users/augustuspendleton/.julia/…
    2 │ 0.191237   lcl|NC_000964.3_cds_NP_387883.1_2  /Users/augustuspendleton/.julia/…
    3 │ 0.0963239  lcl|NC_000964.3_cds_NP_387885.1_4  /Users/augustuspendleton/.julia/…
    4 │ 0.345211   lcl|NC_000964.3_cds_NP_387886.2_5  /Users/augustuspendleton/.julia/…
  ⋮   │     ⋮                      ⋮                                  ⋮
 3799 │ 0.18407    lcl|NC_000964.3_cds_NP_391983.1_…  /Users/augustuspendleton/.julia/…
 3800 │ 0.156161   lcl|NC_000964.3_cds_NP_391984.1_…  /Users/augustuspendleton/.julia/…
 3801 │ 0.235914   lcl|NC_000964.3_cds_NP_391985.1_…  /Users/augustuspendleton/.julia/…
                                                                       3794 rows omitted

```
"""
function scuo(filepath::String, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return scuo(filepath, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function scuo(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Vector{Float64}}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    Threads.@threads for i in 1:len
       @inbounds results[i] = scuo(filepaths[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
    end
    dataframe && return reduce(vcat, results)
    return results
end


function scuo(fasta_seq::String, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # This is our codon count matrix
    names = counts[2] # These are the names of each sequence
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Find amino acid count matrix
    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI) # Calculate normalized frequency of each codon 
    
    Hmax = log10.(dict_deg)

    Oa = map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    Oa = reduce(hcat, Oa)

    Fa = countAA ./ sum(countAA[dict_deg .> 1, :], dims = 1)

    mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    res =  vec(sum(mult, dims = 1))
    dataframe && return DataFrame(SCUO = res, Identifier = names, File = fasta_seq)
    return res
end

"""
    all_cub(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    all_cub(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate all codon usage bias measures at once. Because many measures require the same initial calculations, this is more efficient than calculating them individually.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_seqs`: by default, codon usage bias for each gene is calculated using the whole genome ("self") as a reference subset. If you would like to specify your own subsets to calculate against, such as ribosomal genes, `ref_seqs` takes a named tuple in the form `("subset_name" = Bool[],)`, where `Bool[]` is the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. You can provide multiple reference subsets as separate entries in the named tuple, and `CUBScout` will return the calculated measure using each subset. If providing multiple filepaths and want custom reference sets, `ref_seqs` should be a vector of named tuples corresponding to the vector of filepaths.
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors or named tuples (if the measure uses reference subsets). Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
"""
function all_cub(filepath::String, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    return all_cub(filepath, ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
end

function all_cub(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_seqs = (), rm_start = false, rm_stop = false, threshold = 80, dataframe)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
    end
    if isempty(ref_seqs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = all_cub(filepaths[i], ref_seqs, uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = all_cub(filepaths[i], ref_seqs[i], uniqueI, deg, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end

function all_cub(fasta_seq::String, ref_seqs, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold) # Count codons in each gene 
    count_matrix = counts[1] # This is our codon count matrix
    seqnames = counts[2] # These are our sequence names
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    lengths =  @views transpose(sum(count_matrix, dims = 1)) # Find lengths of each gene (in codons)
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    if isempty(ref_seqs)
        (ref_seqs = (self = fill(true, seqs),)) # If no reference subset provided, make a "self" subset
    else
    ref_seqs = map(x->x[counts[3]], ref_seqs)
    end
    @inbounds countAA = countsbyAA(count_matrix,dict_uniqueI) # Count amino acids in each gene
    @inbounds normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) # Calculate codon frequency within each amino acid within each gene
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs) # Calculate codon frequency within each amino acid across reference subsets
    # Up to this point, all of the measures should have been the same

    # Calculating B
    b_pa = map(x-> x ./ lengths, eachrow(countAA))
    b_pa = transpose(reduce(hcat, b_pa))
    
    b_dts = map(x->abs.((normfreq .- x)), normsetfreqs) # This is good
    b_dts = map(b_dts) do y
        map((x) -> remove_nan(x, 0.0), y)
    end
    b_bas = map(b_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end
    B_result = map(b_bas) do ba
         vec(sum(ba .* b_pa, dims = 1))
    end
    # End calculating B (B is the result)

    # Calculate ENC
    pi_vec =  enc_pi(count_matrix, countAA, seqs, dict_uniqueI)
    enc_fa =  @. (countAA * pi_vec - 1) / (countAA - 1)
    enc_fa[isnan.(enc_fa)] .= 0.0
    ENC_result = (ENC = vec(eFFNc(enc_fa, dict_deg)),)
    # End calculating ENC (ENC is the result)

    # Calculate ENC'
    encp_dts = map(x->(@. (normfreq - x) ^ 2 / x), normsetfreqs) # Calculate Ma for each codon
    encp_dts = map(encp_dts) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    encp_chisums = map(encp_dts) do dt
        chisum = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            chisum[i,:] = row
        end
        chisum
    end
    encp_chisqs = map(x-> x .* countAA, encp_chisums) # This also looks good
    encp_fas = map(encp_chisqs) do chisq
            fa = @. (chisq + countAA - dict_deg) / ((countAA - 1) * dict_deg)
            fa[countAA .< 5] .= 0.0
            fa
        end
    ENCP_result = map(x->vec(eFFNc(x, dict_deg)), encp_fas)
    #End calculating ENC' (ENCP is the result)

    #Start calculating MCB (can use ENCP dts)
    no_counts = count_matrix .<= 0 
    mcb_dts = map(encp_dts) do dt
        dt[no_counts] .= 0
        dt
    end

    mcb_bas = map(mcb_dts) do dt
        ba = Array{Float64}(undef, size(countAA, 1), size(countAA, 2))
        for (i,aa) in enumerate(dict_uniqueI)
            row = @views sum(dt[aa,:],dims = 1)
            ba[i,:] = row
        end
        ba
    end
    A = countAA[dict_deg .> 1,:] .> 0
    MCB_result = map(mcb_bas) do ba
        mat1 = @. ba[dict_deg > 1,:] * log10(countAA[dict_deg > 1, :])
        mat1 = map((x) -> isnan(x) ? 0.0 : x, mat1)
        vec(sum(mat1, dims = 1) ./ sum(A, dims = 1))
    end
        
    # End calculating MCB (MCB is the result)

    # Start calculating MILC
    cor = correction_term(countAA, lengths, dict_deg) # Calculate correction term
    per_codon_mas = map(x->(@. log(normfreq / x) * count_matrix), normsetfreqs) # Calculate Ma for each codon
    fixed_per_codon_mas = map(per_codon_mas) do y
                                map((x) -> remove_nan(x, 0.0), y)
                            end # Remove NaNs
    mas = map(fixed_per_codon_mas) do pcmas
        ma = Array{Float64}(undef, size(countAA, 1), size(countAA, 2)) # Pre-allocate matrix for Ma across each amino acid
        for (i,aa) in enumerate(dict_uniqueI)
            row = 2 .* sum(pcmas[aa,:],dims = 1) # Calculate ma for each amino acid
            ma[i,:] = row
        end
        ma
    end
    MILC_result = map(mas) do ma
        @views vec(([sum(ma, dims = 1)...] ./ lengths) .- cor) # Calculate MILC for each gene
    end
    # End calculating MILC (MILC is the result)

    # Start calculating SCUO
    Ha = scuo_freq(count_matrix, countAA, seqs, dict_uniqueI)
    
    Hmax = log10.(dict_deg)

    Oa = map(x -> (Hmax .- x) ./ Hmax, eachcol(Ha))
    Oa = reduce(hcat, Oa)

    Fa = countAA ./ sum(countAA[dict_deg .> 1, :], dims = 1)

    mult = Oa .* Fa
    mult = map((x) -> isnan(x) ? 0.0 : x, mult)

    SCUO_result = (SCUO = vec(sum(mult, dims = 1)),)

    if dataframe
        milc_df = DataFrame(MILC_result); rename!(milc_df, names(milc_df) .* "_MILC")
        mcb_df = DataFrame(MCB_result); rename!(mcb_df, names(mcb_df) .* "_MCB")
        encp_df = DataFrame(ENCP_result); rename!(encp_df, names(encp_df) .* "_ENCP")
        b_df = DataFrame(B_result); rename!(b_df, names(b_df) .* "_B")
        enc_df = DataFrame(ENC_result)
        scuo_df = DataFrame(SCUO_result) ; scuo_df[!,:Identifier] = seqnames; scuo_df[!,:File] .= fasta_seq
        return hcat(b_df, enc_df, encp_df, mcb_df, milc_df, scuo_df)
    end
    return (B = B_result, ENC = ENC_result, ENCP = ENCP_result, MCB = MCB_result, MILC = MILC_result, SCUO = SCUO_result)

end



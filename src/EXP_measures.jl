
"""
    melp(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    melp(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate MELP from Supek and Vlahovicek, 2005.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `ref_vector`: reference subset, which is required for `melp`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors. Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = melp(example_data_path, ribosomal_genes); # Calculate MELP on example dataset

julia> round.(result.MELP[1:5], digits = 6)
5-element Vector{Float64}:
 0.929414
 1.007671
 0.922357
 0.951239
 1.029531

julia> melp(example_data_path, ribosomal_genes, altstart_codon_dict); # Code TTG and CTG as methionine

julia> melp(example_data_path, ribosomal_genes, rm_start = true); # Remove start codons

julia> melp(example_data_path, ribosomal_genes, dataframe = true); # Get output in dataframe format
```
"""
function melp(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    milcs = milc(filepath, dict, ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector), rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)
    
    if dataframe
    milcs[!,:MELP] = milcs[!,:self] ./ milcs[!,:reference]
    return milcs
    end

    return (MELP = milcs.self ./ milcs.reference, Identifier = milcs.Identifier)
end

function melp(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    ref_tuples = map(x->(self = fill(true, length(x)), reference = x), ref_vectors)
    milcs = milc(filepaths, dict, ref_seqs = ref_tuples, rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)

    if dataframe
    milcs[!,:MELP] = milcs[!,:self] ./ milcs[!,:reference]
    return milcs
    end

    return map(x->(MELP = x.self ./ x.reference, x.Identifier), milcs)
end

"""
    e(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    e(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate E from Karlin and Mrazek, 1996.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `ref_vector`: reference subset, which is required for `e`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors. Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = e(example_data_path, ribosomal_genes); # Calculate E on example dataset

julia> round.(result.E[1:5], digits = 6)
5-element Vector{Float64}:
 0.762317
 1.025839
 0.875954
 0.986498
 1.111275

julia> e(example_data_path, ribosomal_genes, altstart_codon_dict); # Code TTG and CTG as methionine

julia> e(example_data_path, ribosomal_genes, rm_start = true); # Remove start codons

julia> e(example_data_path, ribosomal_genes, dataframe = true); # Get output in dataframe format
```
"""
function e(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    bs = b(filepath, dict, ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector), rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)
    
    if dataframe
    bs[!,:E] = bs[!,:self] ./ bs[!,:reference]
    return bs
    end

    return (E = bs.self ./ bs.reference, Identifier = bs.Identifier)
end

function e(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    ref_tuples = map(x->(self = fill(true, length(x)), reference = x), ref_vectors)
    bs = b(filepaths, dict, ref_seqs = ref_tuples, rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)

    if dataframe
    bs[!,:MELP] = bs[!,:self] ./ bs[!,:reference]
    return bs
    end

    return map(x->(E = x.self ./ x.reference, Identifer = x.Identifier), bs)
end

"""
    cai(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    cai(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate CAI from Sharp and Li, 1987.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `ref_vector`: reference subset, which is required for `cai`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors. Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = cai(example_data_path, ribosomal_genes); # Calculate CAI on example dataset

julia> round.(result.CAI[1:5], digits = 6)
5-element Vector{Float64}:
 0.844967
 0.88548
 0.817348
 1.072675
 0.834179

julia> cai(example_data_path, ribosomal_genes, altstart_codon_dict); # Code TTG and CTG as methionine

julia> cai(example_data_path, ribosomal_genes, rm_start = true); # Remove start codons

julia> cai(example_data_path, ribosomal_genes, dataframe = true); # Get output in dataframe format
```
"""
function cai(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
        aa_names = dict.AA
    end
    return cai(filepath, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
end

function cai(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
        aa_names = dict.AA
    end

    Threads.@threads for i in 1:len
                @inbounds results[i] = cai(filepaths[i], ref_vectors[i], uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
    end
    dataframe && return reduce(vcat, results)
    return results
end


function cai(fasta_seq::String, ref_vector::Vector{Bool}, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, aa_names::Vector{String}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have

    ref_seqs = (self = fill(true, seqs), reference = ref_vector[counts[3]])
    
    countAA = countsbyAA(count_matrix,dict_uniqueI) # This is the same for all measures
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) 
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs)
    max_aa = fill(0.0, length(aa_names))
    map(dict_uniqueI) do aa
        max_aa[aa] .= maximum(normsetfreqs.reference[aa])
    end
    max_aa[max_aa .== 0] .= 0.5
    nodeg =  dict_uniqueI[dict_deg .== 1]
    map(x->count_matrix[x,:] .= 0,nodeg)
    mult = @. log(normfreq / max_aa) * count_matrix
    mult = remove_nan.(mult,0)
    cai_result = vec(exp.(sum(mult, dims = 1) ./ sum(count_matrix, dims = 1)))
    dataframe && return DataFrame(CAI = cai_result, Identifier = names, File = fasta_seq)
    return (CAI = cai_result, Identifier = names)
end


"""
    fop(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    fop(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate FOP from Ikemura, 1981.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `ref_vector`: reference subset, which is required for `fop`. `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors. Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = fop(example_data_path, ribosomal_genes); # Calculate CAI on example dataset

julia> round.(result.FOP[1:5], digits = 6)
5-element Vector{Float64}:
 0.567816
 0.566845
 0.509695
 0.725
 0.653784

julia> fop(example_data_path, ribosomal_genes, altstart_codon_dict); # Code TTG and CTG as methionine

julia> fop(example_data_path, ribosomal_genes, rm_start = true); # Remove start codons

julia> fop(example_data_path, ribosomal_genes, dataframe = true); # Get output in dataframe format
```
"""
function fop(filepath::String, ref_vector::Vector{Bool}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
        aa_names = dict.AA
    end
    return fop(filepath, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
end

function fop(filepaths::Vector{String}, ref_vectors::Vector{Vector{Bool}}, dict::codon_dict = default_codon_dict; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        deg = dict.deg_nostops
        stop_mask = dict.stop_mask
        aa_names = dict.AA_nostops
    else
        uniqueI = dict.uniqueI
        deg = dict.deg
        stop_mask = fill(true,64)
        aa_names = dict.AA
    end
  
    Threads.@threads for i in 1:len
                @inbounds results[i] = fop(filepaths[i], ref_vectors[i], uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
             
    end
    dataframe && return reduce(vcat, results)
    return results
end


function fop(fasta_seq::String, ref_vector::Vector{Bool}, dict_uniqueI::Vector{Vector{Int32}}, dict_deg::Vector{Int32}, stop_mask::Vector{Bool}, aa_names::Vector{String}, rm_start::Bool, threshold::Integer, dataframe::Bool)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have

    ref_seqs = (self = fill(true, seqs), reference = ref_vector[counts[3]])
    
    countAA = countsbyAA(count_matrix,dict_uniqueI) # This is the same for all measures
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) 
    normsetfreqs = map(x->normTotalFreq(count_matrix[:,x], countAA[:,x], dict_uniqueI), ref_seqs)
    max_aa = fill(0.0, length(aa_names))
    map(dict_uniqueI) do aa
        max_aa[aa] .= maximum(normsetfreqs.reference[aa])
    end
    max_aa[max_aa .== 0] .= 0.5
    nodeg =  dict_uniqueI[dict_deg .== 1]
    map(x->count_matrix[x,:] .= 0,nodeg)
    ra = normfreq ./ max_aa
    count2 = copy(count_matrix)
    count2[ra .< 0.9] .= 0
    fops = vec(sum(count2, dims = 1) ./ sum(count_matrix, dims = 1))
    dataframe && return DataFrame(FOP = fops, Identifier = names, File = fasta_seq)
    return (FOP = fops, Identifier = names)

end

"""
    gcb(filepath::String, dict::codon_dict = default_codon_dict; ref_vector = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    gcb(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_vector = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
Calculate GCB from Merkl, 2003.

# Arguments
- `filepath`: path to fasta file of coding sequences (e.g. .fasta, .fna, .fa). There are no quality checks, so it's assumed that each entry is assumed to be an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `filepaths`: vector of paths to fasta files of coding sequences (e.g. .fasta, .fna, .fa). `CUBScout` is multithreaded; if there are multiple threads available, `CUBScout` will allocate a thread for each filepath. As such, providing a vector of paths as an argument will be faster than broadcasting across a vector of paths. Because a single file is only accessed by a single thread, it's never worth using more threads than the total number of files being analyzed. There are no quality checks, so it's assumed that each entry is an individual coding sequence, in the correct frame, without 5' or 3' untranslated regions.
- `dict`: codon dictionary of type `codon_dict`. The standard genetic code is loaded by default, but if necessary you can create your own codon dictionary using `make_codon_dict`
- `ref_vector`: optional reference subset; by default gcb begins calculations using all genes as a seed. If you want to provide a custom reference set, it should be a vector `Bool[]` the same length as the number of sequences in your fasta file, and contains `true` for sequences you want as your reference subset and false for those you don't. You can use `find_seqs()` to generate this vector. If providing multiple filepaths and want custom reference sets, `ref_vectors` should be a vector of vectors corresponding to the vector of filepaths.
- `perc`: percentage of "top hits" which should be used as a reference set in the next iteration. By default set to 0.05. 
- `rm_start`: whether to ignore the first codon of each sequence. Many organisms use alternative start codons such as TTG and CTG, which in other locations would generally code for leucine. There are a few approaches to deal with this. By default, `CUBScout` keeps each start codon and assigns it as though it were any other codon. Of course, this would slightly change leucine's contribution to codon usage bias. If you set `rm_start` to `true`, the first codon of every sequence is simply discarded. This will also affect the gene's length, which means it could be removed if it falls under the threshold. Other CUB packages (such as R's coRdon, alt.init = TRUE), assign all TTG and CTG codons to methionine, regardless of their location. I disagree with this approach from a biological perspective; those codons still code for leucine most of the time they are used. However, if you want matching output as you would get from coRdon, you can supply `altstart_codon_dict` to the `dict` argument, and keep `rm_start` as `false`.
- `rm_stop`: whether to remove stop codons from calculations of codon usage bias.
- `threshold`: minimum length of a gene (in codons) to be used in codon usage bias calculations. By default this is set to 80 codons; any genes less than or equal to that length are discarded. If you want no genes discarded, set `threshold` to 0.
- `dataframe`: whether to format output as a dataframe. By default results are returned as vectors. Setting `dataframe = true` will instead output a dataframe, though at a slight performance cost.

# Examples
```jldoctest
julia> ribosomal_genes = find_seqs(example_data_path, r"ribosomal"); # Get a vector which is true for ribosomal genes

julia> result = gcb(example_data_path); # Calculate GCB on example dataset

julia> round.(result.GCB[1:5], digits = 6)
5-element Vector{Float64}:
 -0.058765
 -0.08659
 -0.005496
 -0.065659
 -0.032062

julia> ribo_result = gcb(example_data_path, ref_vector = ribosomal_genes); # Calculate GCB with ribosomal genes as reference seed example dataset

julia> round.(ribo_result.GCB[1:5], digits = 6)
5-element Vector{Float64}:
 -0.135615
 -0.036687
 -0.169136
 -0.186104
 -0.01653

julia> gcb(example_data_path, altstart_codon_dict); # Code TTG and CTG as methionine

julia> gcb(example_data_path, rm_start = true); # Remove start codons

julia> gcb(example_data_path, dataframe = true); # Get output in dataframe format
```
"""
function gcb(filepath::String, dict::codon_dict = default_codon_dict; ref_vector = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    return gcb(filepath, ref_vector, uniqueI, perc, stop_mask, rm_start, threshold,dataframe)
end

function gcb(filepaths::Vector{String}, dict::codon_dict = default_codon_dict; ref_vectors = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    if isempty(ref_vectors)
        Threads.@threads for i in 1:len
           @inbounds results[i] = gcb(filepaths[i], ref_vectors, uniqueI, perc, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = gcb(filepaths[i], ref_vectors[i], uniqueI, perc, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end


function gcb(fasta_seq::String, refs, dict_uniqueI, perc, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
    counts = count_codons(fasta_seq, rm_start, threshold)# Count codons in each gene 
    count_matrix = counts[1]
    names = counts[2] 
    count_matrix = count_matrix[stop_mask,:] # Remove entries if removing stop codons
    seqs = @views size(count_matrix, 2) # Count how many genes we have
    lengths =  @views transpose(sum(count_matrix, dims = 1)) 

    isempty(refs) ? (seed = fill(true, seqs)) : (seed = refs[counts[3]]) # Make our seed - this will be our initial reference set
    countAA = countsbyAA(count_matrix,dict_uniqueI) # Count 
    normfreq = normFrequency(count_matrix, countAA, seqs, dict_uniqueI) 
    normsetfreq = normTotalFreq(count_matrix[:,seed], countAA[:,seed], dict_uniqueI)
    
    gcb_prev = fill(0.0, seqs)
    iter = 0
    gcb = []
    diff = false
    # Now we'd enter the repeat loop
    while true
        cb = log.(normsetfreq ./ nanmean(normfreq,2))
        cb[normsetfreq .== 0] .= -5
        gcb = vec(vec(sum(count_matrix .* cb, dims = 1)) ./ lengths)
        diff = all(gcb .== gcb_prev)
        if diff | iter > 6
           break
        end
        iter += 1
        gcb_prev = copy(gcb)
        tops = sortperm(gcb, rev = true)[1:convert(Int,trunc(perc*seqs))]
        seed .= false
        seed[tops] .= true
        normsetfreq = normTotalFreq(count_matrix[:,seed], countAA[:,seed], dict_uniqueI)
    end
    dataframe && return DataFrame(GCB = gcb, Identifier = names, File = fasta_seq)
    return (GCB = gcb, Identifier = names)
end


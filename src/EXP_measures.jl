

function melp(filepaths::String, dict::codon_dict, ref_vector::Vector{Bool}; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    milcs = milc(filepaths, dict, ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector), rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)
    
    if dataframe
    milcs[!,:MELP] = milcs[!,:self] ./ milcs[!,:reference]
    return milcs
    end

    return @. milcs.self / milcs.reference
end

function e(filepaths::String, dict::codon_dict, ref_vector::Vector{Bool}; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    bs = b(filepaths, dict, ref_seqs = (self = fill(true, length(ref_vector)), reference = ref_vector), rm_start = rm_start, rm_stop = rm_stop, threshold = threshold, dataframe = dataframe)
    
    if dataframe
    bs[!,:E] = bs[!,:self] ./ bs[!,:reference]
    return bs
    end

    return @. bs.self / bs.reference
end


function cai(filepaths::String, dict::codon_dict, ref_vector::Vector{Bool}; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
    return cai(filepaths, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
end

function cai(filepaths::Vector{String}, dict::codon_dict, ref_vectors; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
    return cai_result
end


# FOP
function fop(filepaths::String, dict::codon_dict, ref_vector::Vector{Bool}; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
    return fop(filepaths, ref_vector, uniqueI, deg, stop_mask, aa_names, rm_start, threshold, dataframe)
end

function fop(filepaths::Vector{String}, dict::codon_dict, ref_vectors; rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
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
    return fops

end

# Finally, GCB
function gcb(filepaths::String, dict::codon_dict; refs = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    return gcb(filepaths, refs, uniqueI, perc, stop_mask, rm_start, threshold,dataframe)
end

function gcb(filepaths::Vector{String}, dict::codon_dict; refs = [], perc = 0.05, rm_start = false, rm_stop = false, threshold = 80, dataframe = false)
    len = length(filepaths)
    results = Vector{Any}(undef, len)
    if rm_stop
        uniqueI = dict.uniqueI_nostops
        stop_mask = dict.stop_mask
    else
        uniqueI = dict.uniqueI
        stop_mask = fill(true,64)
    end
    if isempty(refs)
        Threads.@threads for i in 1:len
           @inbounds results[i] = gcb(filepaths[i], refs, uniqueI, perc, stop_mask, rm_start, threshold, dataframe)
            end
    else
            Threads.@threads for i in 1:len
                @inbounds results[i] = gcb(filepaths[i], refs[i], uniqueI, perc, stop_mask, rm_start, threshold, dataframe)
             end
    end
    dataframe && return reduce(vcat, results)
    return results
end


function gcb(fasta_seq::String, refs, dict_uniqueI::Vector{Vector{Int32}}, perc, stop_mask::Vector{Bool}, rm_start::Bool, threshold::Integer, dataframe)
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
    return gcb
end


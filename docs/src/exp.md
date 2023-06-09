# Expressivity Predictions

## Under default conditions

Expressivity predictions based on codon usage bias can be calculated with the functions `cai()`, `e()`, `fop()`, `gcb()`, and `melp()`. All expressivity functions (besides `gcb`) require two arguments:

- `sequences`: DNA or RNA sequences to be analyzed, which should be coding sequences only. This can take quite a few forms depending on your use case. It can be a path to fasta file of coding sequences (e.g. .fasta, .fna, .fa), or a IO or FASTAReader pointing to these fasta files. It can also be a vector of BioSequences, if you've already brought them into Julia's environment.  If you are analyzing multiple genomes (or sets of sequences), `sequences` could instead be a vector of filepaths, IOStreams, FASTAReaders, or vectors of BioSequences, with each vector corresponding to a genome.
- `ref_vector(s)`: `Vector{Bool}` or `Vector{Vector{Bool}}` identifying reference subsets for each file. Values of `true` should correspond to sequences to be used in the reference subset.

!!! note "Why do expressivity functions accept reference subsets in a different format than codon usage bias functions?"
    You may have noticed that for codon usage bias functions, you need to provide a named tupled of reference sequences, while in expressivity functions, you just need to provide the vector. Why did I make this so complicated for you? Well, for expressivity functions, a reference subset is **required**, and you can only provide a single reference subset. Because this is more strict, the function could be written less flexibly. However, codon usage bias functions can accept no reference subsets, one reference subset, or multiple reference subsets. As such, the named tuple format is necessary to differentiate your input and provide differentiated output.

Let's calculate MELP on our example data, using ribosomal proteins as a reference subset. 

```julia-repl
julia> ribosomal_genes = find_seqs(EXAMPLE_DATA_PATH, r"ribosomal");

julia> melp_result = melp(EXAMPLE_DATA_PATH, ribosomal_genes);

julia> melp_result.MELP
3801-element Vector{Float64}:
 0.9294138732153456
 1.007671319249364
 0.9223573085968517
 0.9512392602630869
 1.0295311265835025
 1.0749743120487463
 ⋮
 1.0147105407479773
 0.9929945778238751
 0.9682178480589456
 0.9651731383865032
 0.8840414848184831
```

The functions `cai`, `e`, `fop`, and `melp` all accept the same arguments. Their optional arguments are the same as codon usage bias functions, including options to specify a custom `CodonDict`, remove start or stop codons, and set a filtering threshold. They also handle multiple files and multi-threading in the same way, and so I do not recommend broadcasting these functions.

## GCB-specific Arguments

Because of the iterative way GCB is calculated, its arguments differ slightly from other expressivity functions. Namely:

- GCB uses a "seed" reference subset. By default, this is set to "self", and so is an optional argument. Custom `ref_vector(s)` can be supplied if so desired, as a keyword argument
- GCB iteratively calculates the GCB measure and then uses the genes with highest GCB values as a reference subset in the next iteration. The `perc` argument specifies what percent of the of genes is used a reference subset. By default, `perc = 0.05`, or 5%. 



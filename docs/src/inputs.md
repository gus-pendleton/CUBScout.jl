# Providing Inputs to CUBScout

## Genomes and CDSs

First, `CUBScout` **only** works with coding sequences. `CUBScout` does not identify ORFs, pause at stop codons, or parse non-nucleotide characters. It is assumed the coding sequences you provide are in-frame and don't contain 5' or 3' untranslated regions. Codons which have non-specific nucleotides, like "W", are skipped.  Sequences with characters outside of those recognized by BioSequences will throw an error. 

Some `CUBScout` functions, like `count_codons`, are meaningful when applied to a single nucleotide sequence. However, most `CUBScout` functions are designed to work at the genome-level, and calculate metrics that rely on comparisons between multiple genes. Specifically, none of the codon usage bias or expressivity functions accept a single nucleotide sequence; all expect to operate across a set of sequences, whether in a fasta file or vector of BioSequences.

## FASTA Files

Most functions in `CUBScout` accept any FASTA-formatted file (e.g. .fa, .fna, .fasta) where each entry corresponds to coding sequences or open readings frames. `CUBScout` accepts either a `String` which is the complete filepath to a fast-formatted file, or objects of type `FASTAReader` or `IO` which point to a fasta-formatted file. There is no significant performance advantage between these three options, unless you already have an `IOStream` or `FASTAReader` open for another purpose. 

## BioSequences

`CUBScout` functions also accept nucleotide sequences from BioSequences (`<:NucSeq`). Keep in mind that most `CUBScout` functions are designed to operate across genomes, and so accept a vector of nucleotide sequences. The vector corresponds to a genome, with each DNA or RNA string corresponding to a coding sequence.

While there is a slight performance advantage in `CUBScout` functions when supplying BioSequences as an input rather than a filepath, supplying filepaths will still be faster than the cumulative time spent reading in a BioSequence and then running a `CUBScout` function. This will also use less memory and so is generally recommended, unless you already have BioSequences loaded into Julia's environment for a separate reason.

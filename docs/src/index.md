```@meta
CurrentModule = CUBScout
```
![header](https://raw.githubusercontent.com/gus-pendleton/CUBScout.jl/main/docs/src/assets/big_pic.png)
# **C**odon **U**sage **B**ias in Julia

[CUBScout](https://github.com/gus-pendleton/CUBScout.jl) helps you work with codons in Julia. You can calculate **C**odon **U**sage **B**ias (CUB) and related expressivity predictions. Currently, `CUBScout` calculates:

- Six measures of codon usage bias:
    - B, from Karlin and Mrazek, 1996
    - ENC, from Wright 1990
    - ENC', from Novembre, 2002
    - MCB, from Urrutia and Hurst, 2001
    - MILC, from Supek and Vlahovicek, 2005
    - SCUO, from Wan et al., 2004
- Five expressivity measures based on codon usage bias:
    - CAI, from Sharp and Li, 1987
    - E, from Karlin and Mrazek, 1996
    - FOP, from Ikemura, 1981
    - GCB, from Merkl, 2003
    - MELP, from Supek and Vlahovicek, 2005

`CUBScout` is based off of the fabulous [coRdon](https://www.bioconductor.org/packages/release/bioc/html/coRdon.html) package in R by Anamaria Elek, Maja Kuzman, and Kristian Vlahovicek. I am grateful for their clear code and would encourage you to cite coRdon as well when using `CUBScout`.



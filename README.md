# ITensorPySCF

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://mtfishman.github.io/ITensorPySCF.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://mtfishman.github.io/ITensorPySCF.jl/dev)
[![Build Status](https://github.com/mtfishman/ITensorPySCF.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/mtfishman/ITensorPySCF.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/mtfishman/ITensorPySCF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mtfishman/ITensorPySCF.jl)
[![Code Style: Blue](https://img.shields.io/badge/code%20style-blue-4495d1.svg)](https://github.com/invenia/BlueStyle)

A package for interface [ITensors.jl](https://github.com/ITensor/ITensors.jl) and [PySCF](https://github.com/pyscf/pyscf), using PyCall.

The goal is to use PySCF to provide one and two body integrals and mean-field states to run quantum chemistry tensor network calculations like DMRG.

Install with:
```julia
julia> ] add https://github.com/ITensor/ITensorGaussianMPS.jl#main

julia> ] add https://github.com/mtfishman/ITensorPySCF.jl#main
```

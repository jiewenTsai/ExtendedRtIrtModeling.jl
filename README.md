# ExtendedRtIrtModeling.jl

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://jiewenTsai.github.io/ExtendedRtIrtModeling.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://jiewenTsai.github.io/ExtendedRtIrtModeling.jl/dev)
[![Build Status](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/workflows/Test/badge.svg)](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions)
[![Test workflow status](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/actions/workflows/Docs.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/jiewenTsai/ExtendedRtIrtModeling.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/jiewenTsai/ExtendedRtIrtModeling.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](CODE_OF_CONDUCT.md)
[![All Contributors](https://img.shields.io/github/all-contributors/jiewenTsai/ExtendedRtIrtModeling.jl?labelColor=5e1ec7&color=c0ffee&style=flat-square)](#contributors)
[![BestieTemplate](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/JuliaBesties/BestieTemplate.jl/main/docs/src/assets/badge.json)](https://github.com/JuliaBesties/BestieTemplate.jl)

## Overview

Four main RtIrt models are provided in this package,

- `GibbsMlIrt`
- `GibbsRtIrt`
- `GibbsRtIrtQuantile`

These three models default to account for covariate variables (e.g., latent regression, and latent structure). If you need only a measurement model, you can use the null model.

- `GibbsRtIrtNull`


## Installation

You can download `ExtendedRtIrtModeling` directly from julia.

```julia
using Pkg
Pkg.add("ExtendedRtIrtModeling")
```

or 

```julia
] add ExtendedRtIrtModeling
```


This package isn’t registered in Julia yet, so you must download it from GitHub.

```julia
using Pkg
Pkg.add(url="https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl")
```

or 

```julia
] add "https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl"
```

## Get Started

Here is a simulation study example.

```julia
using ExtendedRtIrtModeling

## creat a toy data
Cond = setCond(nSubj=1000, nItem=15)
truePara = setTrueParaMlIrt(Cond)
Data = setDataMlIrt(Cond, truePara)

## build a model and sample it!
MCMC = GibbsMlIrt(Cond, Data=Data, truePara=truePara)
sample!(MCMC)

## check the parameter recovery
getRmse(MCMC.truePara.b, MCMC.Post.mean.b)

```

If you have a data set to analyze, you can follow the following way,

```julia
using ExtendedRtIrtModeling
using CSV, DataFrames

## import your data set
yourData = CSV.read("yourData.csv", DataFrame)
Cond = setCond(qRa=0.85, qRt=0.85, nChain=3, nIter=3000)
Data = InputData(
    Y=Matrix(yourData[:,1:15]),
    T=exp.(Matrix(yourData[:,16:30])),
    X=Matrix(yourData[:,31:33])
)

## build a model and sample it!
MCMC = GibbsRtIrtQuantile(Cond, Data=Data)
sample!(MCMC)


MCMC.Post.mean.Σp
MCMC.Post.mean.β

```



## How to Cite

If you use ExtendedRtIrtModeling.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/jiewenTsai/ExtendedRtIrtModeling.jl/blob/main/CITATION.cff).


## Contributing

If you want to make contributions of any kind, please first that a look into our [contributing guide directly on GitHub](docs/src/90-contributing.md) or the [contributing page on the website](https://jiewenTsai.github.io/ExtendedRtIrtModeling.jl/dev/90-contributing/).


---

### Contributors

<!-- ALL-CONTRIBUTORS-LIST:START - Do not remove or modify this section -->
<!-- prettier-ignore-start -->
<!-- markdownlint-disable -->

<!-- markdownlint-restore -->
<!-- prettier-ignore-end -->

<!-- ALL-CONTRIBUTORS-LIST:END -->


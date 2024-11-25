
include("src/ExtendedRtIrtModeling.jl")
import .ExtendedRtIrtModeling as ex

using Distributions,
    LinearAlgebra,
    Plots,
    HTTP,
    CSV,
    DataFrames


Cond = ex.setCond(nChain=2, nIter=3000, nSubj=1000, nItem=10)
truePara = ex.setTrueParaRtIrt(Cond, trueCorr=0.5)
Data = ex.setDataRtIrtNull(Cond, truePara)
#Data2 = ex.setDataMlIrt(Cond, truePara)


MCMC = ex.GibbsRtIrtCross(Cond, Data=Data, truePara=truePara)
ex.sample!(MCMC)

ex.coef(MCMC)


ex.getRmse(MCMC.truePara.a, MCMC.Post.mean.a)
ex.getBias(MCMC.truePara.a, MCMC.Post.mean.a)
[MCMC.truePara.a MCMC.Post.mean.a]
ex.getRmse(MCMC.truePara.b, MCMC.Post.mean.b)
ex.getBias(MCMC.truePara.b, MCMC.Post.mean.b)
[MCMC.truePara.b MCMC.Post.mean.b]
ex.getRmse(MCMC.truePara.λ, MCMC.Post.mean.λ)
ex.getRmse(MCMC.truePara.σ²t, MCMC.Post.mean.σ²t)
ex.getRmse(vec(MCMC.truePara.β), vec(MCMC.Post.mean.β))



histogram(MCMC.truePara.θ)
histogram!(MCMC.Post.mean.θ ./ 0.85)
(MCMC.Post.mean.θ) |> mean
MCMC.Post.mean.b |> mean


(sum(MCMC.truePara.a - MCMC.Post.mean.a) / 10)
sqrt(sum((MCMC.truePara.a .- MCMC.Post.mean.a).^2) / 10)

sqrt(sum((MCMC.truePara.a .- MCMC.Post.mean.a).^2) / 10)

length(MCMC.truePara.a)

# =================
# Real data
# =================

demoHttp = ("https://raw.githubusercontent.com/jiewenTsai/ExtendedRtIrtModeling.jl/refs/heads/main/data/demo.csv")
Demo = CSV.read(HTTP.get(demoHttp).body, DataFrame)
Data7 = ex.InputData(
    Y=Matrix(Demo[:,2:11]),
    T=Matrix(exp.(Demo[:,12:21])),
    X=Matrix(Demo[:,22:25])
)


#begin
Cond7_50 = ex.setCond(
    nSubj=300, nItem=10, nFeat=4, nChain=3, nIter=1500,
    nThin=3,
    qRt=0.5,
)
MCMC7 = ex.GibbsRtIrtNull(Cond7_50, Data=Data7)
ex.sample!(MCMC7, cov2one=true)
#end

ex.coef(MCMC7)
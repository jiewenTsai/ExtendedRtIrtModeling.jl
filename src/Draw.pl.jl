### A Pluto.jl notebook ###
# v0.20.3

#using Markdown
#using InteractiveUtils

# ╔═╡ f134213b-7eeb-4216-bfdc-1ced3c3c3918
using PlutoUI

# ╔═╡ ab165b5e-9f8a-11ef-231d-6f8f057f3cf8
md"""

# *Draw*-functions

The *draw-*functions we provided are as follows,

- 

"""

# ╔═╡ fb72b3fe-c677-4969-b262-6a63fede0752
TableOfContents()

# ╔═╡ 307e75db-26f9-48e9-b7c1-1ac462dc3e9d
md"""
## IRT model
"""

# ╔═╡ 85b358b6-d32a-4a6e-b61e-bbba8d02fc46
md"### Random"

# ╔═╡ 9d8d02c3-4fca-4f22-9a5b-d32cdc632cac
"""
	drawRaPgRandomVariable(Para) --> ω
"""
function drawRaPgRandomVariable(Para)
	η = Para.a' .* (Para.θ .- Para.b')
	ω = rand.(PolyaGammaPSWSampler.(1, η))
    return ω
end

# ╔═╡ 5388c6fe-29c4-421c-82b4-c1bb9468683f
md"### Subject parameter - Accuracy"

# ╔═╡ 5fd59c14-e0c8-4c22-86d9-5e4b3e7e48f1
"""
    drawSubjAbility(Cond,Data,Para) --> θ  
"""
function drawSubjAbility(Cond,Data,Para)
    x = [ones(Cond.nSubj) Data.X]

	θμ₀ = x * Para.β[:,1] 
	θσ₀² = 1. #Para.Σp[1,1]

    parV = 1 ./(1 ./θσ₀² .+ sum(Para.a'.^2 .* Para.ω, dims=2))
	parM = parV .* (θμ₀ ./θσ₀² .+ sum( Para.a' .* (Data.κ .+ Para.a' .* Para.b' .* Para.ω ), dims=2 ))
	θ = rand.(Normal.(parM, sqrt.(parV))) 

	#prodA = prod(Para.a)^(1/Cond.nItem)
	#θNew = θ .* prodA
	return θ
end

# ╔═╡ 67ac71e7-9250-4130-95f7-351e62efd925
"""
"""
function drawSubjAbilityNull(Cond,Data,Para)
    #x = [ones(Cond.nSubj) Data.X]
	θμ₀ = 0.
	θσ₀² = 1. #Para.Σp[1,1] 

    parV = 1 ./(1 ./θσ₀² .+ sum(Para.a'.^2 .* Para.ω, dims=2))
	parM = parV .* (θμ₀ ./θσ₀² .+ sum( Para.a' .* (Data.κ .+ Para.a' .* Para.b' .* Para.ω ), dims=2 ))
	θ = rand.(Normal.(parM, sqrt.(parV))) 

	return θ
end

# ╔═╡ de569bf7-a393-48a3-89a1-e7b1bec148ff
md"### Item parameters"

# ╔═╡ 4a977b7d-98e0-4a24-9889-db6bafcba040
"""
"""
function drawItemDiscrimination(Data, Para; μa₀=1., σa₀=1e+10)
	parV = 1 ./ (1/σa₀.^2 .+ sum( (Para.θ .- Para.b').^2 .*  Para.ω, dims=1))
    parM = parV .* ( μa₀/σa₀.^2  .+  sum( Data.κ .* (Para.θ .- Para.b'), dims=1) )
    a = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf)  ) 
    return a'
end

# ╔═╡ 33582955-b9e0-406f-b0ad-995b2d630236
"""
"""
function drawItemDifficulty(Data, Para; μb₀=0., σb₀=1e+10)
	parV = 1 ./ (1/σb₀^2 .+  sum( Para.a'.^2 .* Para.ω, dims=1)')
    parM = parV .* ( μb₀/σb₀^2 .- sum( Para.a' .* ( Data.κ .- Para.θ*Para.a'.*Para.ω), dims=1) )'
    b = rand.(Normal.(parM, sqrt.(parV))) 
    return b
end


# ╔═╡ 63063231-7e9a-4d13-b294-1f4f2939e867
md"""
## RT model
"""

# ╔═╡ 7a2c7239-a78f-4632-9905-959fb79ea50b
md"### Subject parameter - Speed"

# ╔═╡ 6e564b7e-313b-405d-b04c-b41bd15e7fc2
"""
"""
function drawSubjSpeedNull(Cond,Data,Para )
	ζμ₀ = 0. #x * Para.β[:,2] 
	ζσ₀² = 1. #Para.Σp[2,2] 
	
    parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./  Para.σ²t', dims=2) )
    parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum( (Para.λ' .- Data.logT) ./ Para.σ²t', dims=2))  
    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 638123de-6b80-4a36-8b0e-f3009daa9bbe
"""
"""
function drawSubjSpeed(Cond,Data,Para )
	x = [ones(Cond.nSubj) Data.X]
	ζμ₀ = x * Para.β[:,2] 
	ζσ₀² = Para.Σp[2,2]  

    parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./  Para.σ²t', dims=2) )
    parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum( (Para.λ' .- Data.logT) ./ Para.σ²t', dims=2))  
    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 334cce4b-4cfc-475c-9c92-5ecd7ec3c6b1
"""
	drawSubjSpeedLatent(Cond,Data,Para) --> ζ
"""
function drawSubjSpeedLatent(Cond,Data,Para )
	x = [ones(Cond.nSubj) Data.X Para.θ]
	ζμ₀ = x * Para.β 
	ζσ₀² = Para.Σp[2,2]  

    parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./  Para.σ²t', dims=2) )
    parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum( (Para.λ' .- Data.logT) ./ Para.σ²t', dims=2))  
    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 13fd37a1-52cb-4938-a4e4-6eba043824df
"""
"""
function drawSubjSpeedLatentQr(Cond,Data,Para )
	x = [ones(Cond.nSubj) Data.X Para.θ]
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))

	## think more about this part
	ζμ₀ =  x * Para.β + k1Rt * Para.ν
	ζσ₀² = Para.Σp[2,2] .* (k2Rt * Para.ν)
	
    parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./  Para.σ²t', dims=2) )
    parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum( (Para.λ' .- Data.logT) ./ Para.σ²t', dims=2))  
    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 3e6213f4-d6fc-40de-a71d-053740dd068b
"""
"""
function drawSubjSpeedCross(Cond,Data,Para )
	ζμ₀ = 0.
	ζσ₀² = Para.Σp[2,2] 

	parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./ Para.σ²t', dims=2))
    parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum((Para.λ' .- Data.logT .- Para.θ * Para.ρ') ./ Para.σ²t', dims=2))
    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 1aee48af-2dbf-4100-80c0-d4cf7bb5c02c
"""
"""
function drawSubjSpeedCrossQr(Cond,Data,Para )
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν

	ζμ₀ = 0.
	ζσ₀² = Para.Σp[2,2] 
	
    parV = 1 ./ (1 ./ ζσ₀² .+ sum( 1. ./ (Para.σ²t' .* k2e ), dims=2))
	parM = parV .* (ζμ₀ ./ ζσ₀² .+ sum( (Para.λ' .- Data.logT .- Para.θ * Para.ρ' .+ k1e) ./ (Para.σ²t' .* k2e ), dims=2))

    ζ = rand.(Normal.(parM, sqrt.(parV)))
    return ζ
end

# ╔═╡ 1ed35050-a370-482e-8aad-ab025a6c30e4
md"### Item parameters"

# ╔═╡ 0b10b047-034d-4f5d-9832-a058fb3cd8fa
"""
    drawItemIntensity(Cond,Data,Para; μλ=mean(Data.logT), σλ=1e+10) --> λ
"""
function drawItemIntensity(Cond,Data,Para; μλ=mean(Data.logT), σλ=std(Data.logT))
    parV = 1 ./(1/σλ^2 .+ Cond.nSubj ./ Para.σ²t)
    parM = parV .* (μλ/σλ^2 .+ sum(Data.logT .+ Para.ζ, dims=1) ./ Para.σ²t')'
    λ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf))
    return λ
end

# ╔═╡ 002c3afe-12a8-480b-8bb0-b974a5dfcff4
"""
"""
function drawItemIntensityCross(Cond,Data,Para; μλ=mean(Data.logT), σλ=std(Data.logT))
	parV = 1 ./( 1/σλ^2 .+ sum(Cond.nSubj ./ Para.σ²t', dims=1))
    parM = parV .* (μλ/σλ^2 .+ sum( ( (Data.logT .+ Para.ζ .+ Para.θ * Para.ρ') ./ Para.σ²t'), dims=1))
    λ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf))

    return λ'
end

# ╔═╡ 923f8ff9-d0bb-4e4f-932d-0da8600b9318
"""
	drawItemIntensityCrossQr

Note. 
"""
function drawItemIntensityCrossQr(Cond,Data,Para; μλ=mean(Data.logT), σλ=std(Data.logT))

	k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν
	
	parV = 1 ./( 1/σλ^2 .+ sum( 1. ./ (Para.σ²t' .* k2e), dims=1 ))
	parM = parV .* (μλ/σλ^2 .+ ( sum( (Data.logT .+ Para.ζ .+ Para.θ * Para.ρ' .- k1e) ./ (Para.σ²t' .* k2e), dims=1)))
    λ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf))

    return λ'
end

# ╔═╡ bdf1e8e9-6c01-4760-be17-a91a24832700
"""
    drawItemTimeResidual(Cond,Data,Para;δa=1e-10, δb=1e-10) --> σ²t
"""
function drawItemTimeResidual(Cond,Data,Para;δa=1e-10, δb=1e-10)
    parA = δa + Cond.nSubj ./2
    parB = δb .+ sum( (Data.logT .- Para.λ' .+ Para.ζ).^2, dims=1) ./2
    σ²t = rand.(InverseGamma.(parA, parB'))
    return σ²t
end

# ╔═╡ 84e5e58b-02cc-4fd0-8a60-67b393f8b8fe
"""
"""
function drawItemTimeResidualCross(Cond,Data,Para;δa=1e-10, δb=1e-10)
    parA = δa + Cond.nSubj/2
    parB = δb .+ sum( (Data.logT .- Para.λ' .+ Para.ζ .+ Para.θ * Para.ρ').^2, dims=1) ./2
    σ²t = rand.(InverseGamma.(parA, parB'))

    return σ²t
end

# ╔═╡ 982acb19-e7b1-4cd6-a98d-940d1a4c1565
"""
"""
function drawItemTimeResidualCrossQr(Cond,Data,Para;δa=1e-10, δb=1e-10)
	k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν

    parA = δa + Cond.nSubj*3 /2 
	parB = δb .+ sum( (Data.logT .- Para.λ' .+ Para.ζ .+ Para.θ * Para.ρ' .- k1e).^2 ./ (2 .* k2e), dims=1) .+ sum(Para.ν, dims=1)
    σ²t = rand.(InverseGamma.(parA, parB'))
    return σ²t
end

# ╔═╡ 5293be72-3fb3-47b1-9d9d-b5f6228605f5
md"""

## Structural model

"""

# ╔═╡ 0bd6c7df-2c47-4295-84eb-71fa21306e60
md"### Weights"

# ╔═╡ 4e053773-b880-4a19-a5c2-9f86c74485ca
"""
"""
function drawQrWeightsCrossQr(Cond,Data,Para)
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
	
	parA =  abs.(Data.logT .- Para.λ' .+ Para.ζ .+ Para.θ * Para.ρ') ./ sqrt.(Para.σ²t' .* k2Rt )
	parB =  sqrt.(2 * k2Rt .+ k1Rt.^2) ./ sqrt.(Para.σ²t' .* k2Rt)

	
    ν =  1 ./ rand.(InverseGaussian.( (parB ./ parA), parB.^2  ))
	#ν = 1 ./ rand.(InverseGaussian.( parMu, parLam  ))
	#ν = rand.(ge.GeneralizedInverseGaussian.(0.5,parB,parA ))
	#ν = 1 ./ rand.(ge.GeneralizedInverseGaussian.(-0.5,parA,parB ))


    return  ν
end

# ╔═╡ 39381143-5491-4ed8-9874-805c00050740
"""
"""
function drawQrWeightsLatentQr(Cond,Data,Para)
	x = [ones(Cond.nSubj) Data.X Para.θ]
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))

	parA =  abs.(Para.ζ .- x * Para.β) ./ sqrt.(Para.Σp[2,2] .* k2Rt)
	parB =  sqrt.(2 * k2Rt .+ k1Rt.^2) ./ sqrt.(Para.Σp[2,2] .* k2Rt) 

	ν =  1 ./ rand.(InverseGaussian.( (parB ./ parA), parB.^2  ))
    #ν =  1 ./ rand.(InverseGaussian.(parMu  , parLam ))
	#ν = rand.(ge.GeneralizedInverseGaussian.(0.5,parB,parA))
	#ν = 1 ./ rand.(ge.GeneralizedInverseGaussian.(-0.5,parA,parB ))

    return ν
end

# ╔═╡ a14d81ba-8d73-446c-b3be-6691da980043
md"### β Coefficients"

# ╔═╡ 66562193-7978-4fb8-b8bd-3430ec74d30a
"""
"""
function getSubjCoefficientsMlIrt(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )
	η = Para.θ
	x = [ones(Cond.nSubj) Data.X]

	β = x'x \ x'*η
    return β
end

# ╔═╡ 617971fa-2469-4155-9ead-1ccb78c372f2
"""
"""
function getSubjCoefficients(Cond, Data, Para )
	#prodA = (prod(Para.a)^(1/Cond.nItem))

	η = [Para.θ Para.ζ]
	x = [ones(Cond.nSubj) Data.X]
	invΩ = inv(Para.Σp)
	#β = x'x \ x'*η
	β = (invΩ ⊗ x'x) \ vec(x'*η * invΩ')
	#β = (invΩ ⊗ x'x) \ vec(x'*η )
	#β = Data.X'Data.X \ (Data.X' * η * inv(Para.Σp))

	β = reshape(β, Cond.nFeat+1,2)
	return β
end

# ╔═╡ 62e6e5bc-ba4f-4c81-94d5-3514b4cc0f5c
"""
"""
function drawSubjCoefficients(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )

    η = [Para.θ Para.ζ]
    x = [ones(Cond.nSubj) Data.X]

    invΩ = inv(Symmetric(Para.Σp))
	parV = inv( 1/σβ₀^2 .+  invΩ ⊗ x'x) 
    parM = parV * ( μβ₀/σβ₀^2 .+ vec(x'* η * invΩ') )

    ## Generate random coefficients
    β = parM .+ cholesky(Symmetric(parV)).L * randn((Cond.nFeat+1)*2)

    return reshape(β, Cond.nFeat+1,2)
end

# ╔═╡ a831605b-f376-491b-adff-388a7ac553df
"""
    drawSubjCoefficientsLatentQr(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 ) --> β
"""
function drawSubjCoefficientsLatent(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )
    x = [ones(Cond.nSubj) Data.X Para.θ]
    """
    parV = 1 ./ (1/σβ₀^2 .+ sum(x.^2 ./ Para.Σp[2,2], dims=1))
    parM = parV .* (μβ₀/σβ₀^2 .+ sum(x .* Para.ζ ./ Para.Σp[2,2], dims=1))   
    β = rand.(Normal.(vec(parM), sqrt.(vec(parV))))
    """

	#β = x'x \ x'Para.ζ

	invΩ = inv(Para.Σp[2,2])
	parV = inv( 1/σβ₀^2 .+  invΩ * x'x) 
    parM = parV * ( μβ₀/σβ₀^2 .+ vec(x'* Para.ζ * invΩ') )

	β = parM .+ cholesky(Symmetric(parV)).L * randn((Cond.nFeat+2))

    return β
end

# ╔═╡ e4c42f0d-a7a2-4877-bbf4-15f8fdfb310e
md"""
!!! info
	這個 drawSubjCoefficientsLatentQr 的問題還沒解決，目前只是先用 get- 擋著。
"""

# ╔═╡ 44923bdc-3ac9-4e3b-887d-b1b771c2834b
"""
    drawSubjCoefficientsLatentQr(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 ) --> β
"""
function drawSubjCoefficientsLatentQr(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )
    x = [ones(Cond.nSubj) Data.X Para.θ]
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν

	parV = 1. ./ (1/σβ₀^2 .+ sum(x.^2  ./ (k2e .* Para.Σp[2,2] ), dims=1))
    parM = parV .* (μβ₀/σβ₀^2 .+ sum(x .* (Para.ζ .- k1e) ./ (k2e .* Para.Σp[2,2] ), dims=1))   
    β = rand.(Normal.(vec(parM), sqrt.(vec(parV))))

    return β
end

# ╔═╡ 2246e0f5-b367-4818-b06c-c921d84d354c
"""
    getSubjCoefficientsLatentQr(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 ) --> β
"""
function getSubjCoefficientsLatentQr(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )
    x = [ones(Cond.nSubj) Data.X Para.θ]
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν
	k2eΣp = 1 ./ (Para.Σp[2,2] .* k2e)
	#k2eΣp =  inv.(Para.Σp[2,2])

	β = (k2eΣp ⊗ x'x)  \ (vec(x' * (Para.ζ .- k1e) * k2eΣp'))
	
    return β
end

# ╔═╡ 30dd078c-b514-4130-aec8-5c04f63c4b02
"""
"""
function drawSubjCorrCross(Cond,Data,Para;μρ=0., σρ=1e+10)
	parV = 1 ./( 1/σρ^2 .+ sum(Para.θ.^2 ./ Para.σ²t', dims=1 ))
	parM = parV .* (μρ/σρ^2 .+  sum( (Para.θ .* (Para.λ' .- Para.ζ .- Data.logT)) ./ Para.σ²t' , dims=1))
	ρ = rand.(Normal.(parM, sqrt.(parV)))

    return ρ'
end

# ╔═╡ 1728fba4-1bfb-4630-a5da-41ff2df84f03
"""
"""
function drawSubjCorrCrossQr(Cond,Data,Para;μρ=0., σρ=1e+10)
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν

	parV = 1 ./( 1/σρ^2 .+ sum(Para.θ.^2 ./ (Para.σ²t' .* k2e ), dims=1 ))
	parM = parV .* (μρ/σρ^2 .+  sum( (Para.θ .* (Para.λ' .- Para.ζ .- Data.logT .+ k1e)) ./ (Para.σ²t' .* k2e) , dims=1))
	ρ = rand.(Normal.(parM, sqrt.(parV)))

    return ρ'
end

# ╔═╡ d3d8b090-e9d8-440d-adfd-222a5226a37f
md"### Σ Covariance"

# ╔═╡ 902dc1c9-e1ff-4044-a01c-76141dd457d6
"""
    drawSubjCovariance(;θ, τ, nSubj)
    (New)
"""
function drawSubjCovariance(Cond, Data, Para, cov2one)
    η = [Para.θ Para.ζ]
    x = [ones(Cond.nSubj) Data.X]
	e = η .- x * Para.β 

    ee = e'e
    s =  rand(InverseWishart(Cond.nSubj+3, ee .+ I(2)) )

	if cov2one == true
		s = diagm([s[1,1].^-0.5,1.]) * s * diagm([s[1,1].^-0.5, 1.])
		s = diagm([1., s[2,2].^-0.5]) * s * diagm([1., s[2,2].^-0.5])
		s[1,1] = s[2,2] =  1.
	end


    return s
end

# ╔═╡ 346430bf-7397-426b-84a9-b322c469c89d
"""
    drawSubjCovariance(;θ, τ, nSubj)
    (New)
"""
function drawSubjCovarianceNull(Cond, Data, Para, cov2one)
    η = [Para.θ Para.ζ]
    ee = η'η
	#s = e'e ./ Cond.nSubj
	s =  rand(InverseWishart(Cond.nSubj+3, ee + I(2)))

	if cov2one == true
		s = diagm([s[1,1].^-0.5,1.]) * s * diagm([s[1,1].^-0.5, 1.])
		s = diagm([1., s[2,2].^-0.5]) * s * diagm([1., s[2,2].^-0.5])
		s[1,1] = s[2,2] =  1.
	end

    return s #Σp
end

# ╔═╡ 31b22562-2177-4f23-9320-ce88686b0ddc
"""
    drawSubjCovariance(;θ, τ, nSubj)
    (New)
"""
function drawSubjCovarianceCross(Cond, Data, Para, cov2one;  δa = 1e-10, δb = 1e-10)
    parA = δa + Cond.nSubj/2
    #parB = δb .+ sum(quantileCheck.(Para.ζ .- x * Para.β, Cond.qRt))/2
    parB = δb .+ sum((Para.ζ).^2)/2
    s =  rand(InverseGamma(parA, parB))
    Σp = [1. 0.; 0. s]

	if cov2one == true
		Σp = diagm([Σp[1,1].^-0.5,1.]) * Σp * diagm([Σp[1,1].^-0.5, 1.])
		Σp = diagm([1., Σp[2,2].^-0.5]) * Σp * diagm([1., Σp[2,2].^-0.5])
		Σp[1,1] = Σp[2,2] =  1.
	end


    return Σp
end

# ╔═╡ a8e5c5df-8def-4498-95f8-2deab8b15366
"""
    drawSubjCovarianceLatentQr(Cond, Data, Para) --> Σp
"""
function drawSubjCovarianceLatent(Cond, Data, Para, cov2one; δa = 1e-10, δb = 1e-10)
    x = [ones(Cond.nSubj) Data.X Para.θ]

    parA = δa + Cond.nSubj/2
    #parB = δb .+ sum(quantileCheck.(Para.ζ .- x * Para.β, Cond.qRt))/2
    parB = δb .+ sum((Para.ζ .- x * Para.β).^2)/2
    s =  rand(InverseGamma(parA, parB))
    Σp = [1. 0.; 0. s]

	if cov2one == true
		Σp = diagm([Σp[1,1].^-0.5,1.]) * Σp * diagm([Σp[1,1].^-0.5, 1.])
		Σp = diagm([1., Σp[2,2].^-0.5]) * Σp * diagm([1., Σp[2,2].^-0.5])
		Σp[1,1] = Σp[2,2] =  1.
	end

    return Σp
end

# ╔═╡ ac6dd3a5-9b1a-40d8-bb9d-0582975dc319
"""
    drawSubjCovarianceLatentQr(Cond, Data, Para) --> Σp
"""
function drawSubjCovarianceLatentQr(Cond, Data, Para, cov2one; δa = 1e-10, δb = 1e-10)
    x = [ones(Cond.nSubj) Data.X Para.θ]
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = k1Rt.*Para.ν
	k2e = k2Rt.*Para.ν

	
    parA = δa + Cond.nSubj*3 /2 
	parB = δb .+ sum((Para.ζ .- x * Para.β .- k1Rt * Para.ν).^2  / (2*k2e)) .+ sum(Para.ν)

	s =  rand(InverseGamma(parA, parB))
    Σp = [1. 0.; 0. s]

	if cov2one == true
		Σp = diagm([Σp[1,1].^-0.5,1.]) * Σp * diagm([Σp[1,1].^-0.5, 1.])
		Σp = diagm([1., Σp[2,2].^-0.5]) * Σp * diagm([1., Σp[2,2].^-0.5])
		Σp[1,1] = Σp[2,2] =  1.
	end

    return Σp
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
PlutoUI = "~0.7.60"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "c1674f662899f5bfc062df83020732df21a649e9"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "b6d6bfdd7ce25b0f9b2f6b3dd56b2673a66c8770"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.5"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "eba4810d5e6a01f612b948c9fa94f905b49087b0"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.60"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.Tricks]]
git-tree-sha1 = "7822b97e99a1672bfb1b49b668a6d46d58d8cbcb"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.9"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─ab165b5e-9f8a-11ef-231d-6f8f057f3cf8
# ╠═f134213b-7eeb-4216-bfdc-1ced3c3c3918
# ╠═fb72b3fe-c677-4969-b262-6a63fede0752
# ╟─307e75db-26f9-48e9-b7c1-1ac462dc3e9d
# ╟─85b358b6-d32a-4a6e-b61e-bbba8d02fc46
# ╠═9d8d02c3-4fca-4f22-9a5b-d32cdc632cac
# ╟─5388c6fe-29c4-421c-82b4-c1bb9468683f
# ╠═5fd59c14-e0c8-4c22-86d9-5e4b3e7e48f1
# ╠═67ac71e7-9250-4130-95f7-351e62efd925
# ╟─de569bf7-a393-48a3-89a1-e7b1bec148ff
# ╟─4a977b7d-98e0-4a24-9889-db6bafcba040
# ╟─33582955-b9e0-406f-b0ad-995b2d630236
# ╟─63063231-7e9a-4d13-b294-1f4f2939e867
# ╟─7a2c7239-a78f-4632-9905-959fb79ea50b
# ╠═6e564b7e-313b-405d-b04c-b41bd15e7fc2
# ╠═638123de-6b80-4a36-8b0e-f3009daa9bbe
# ╠═334cce4b-4cfc-475c-9c92-5ecd7ec3c6b1
# ╠═13fd37a1-52cb-4938-a4e4-6eba043824df
# ╠═3e6213f4-d6fc-40de-a71d-053740dd068b
# ╠═1aee48af-2dbf-4100-80c0-d4cf7bb5c02c
# ╟─1ed35050-a370-482e-8aad-ab025a6c30e4
# ╟─0b10b047-034d-4f5d-9832-a058fb3cd8fa
# ╟─002c3afe-12a8-480b-8bb0-b974a5dfcff4
# ╠═923f8ff9-d0bb-4e4f-932d-0da8600b9318
# ╟─bdf1e8e9-6c01-4760-be17-a91a24832700
# ╟─84e5e58b-02cc-4fd0-8a60-67b393f8b8fe
# ╠═982acb19-e7b1-4cd6-a98d-940d1a4c1565
# ╟─5293be72-3fb3-47b1-9d9d-b5f6228605f5
# ╟─0bd6c7df-2c47-4295-84eb-71fa21306e60
# ╠═4e053773-b880-4a19-a5c2-9f86c74485ca
# ╠═39381143-5491-4ed8-9874-805c00050740
# ╟─a14d81ba-8d73-446c-b3be-6691da980043
# ╠═66562193-7978-4fb8-b8bd-3430ec74d30a
# ╠═617971fa-2469-4155-9ead-1ccb78c372f2
# ╠═62e6e5bc-ba4f-4c81-94d5-3514b4cc0f5c
# ╠═a831605b-f376-491b-adff-388a7ac553df
# ╠═e4c42f0d-a7a2-4877-bbf4-15f8fdfb310e
# ╠═44923bdc-3ac9-4e3b-887d-b1b771c2834b
# ╠═2246e0f5-b367-4818-b06c-c921d84d354c
# ╟─30dd078c-b514-4130-aec8-5c04f63c4b02
# ╟─1728fba4-1bfb-4630-a5da-41ff2df84f03
# ╟─d3d8b090-e9d8-440d-adfd-222a5226a37f
# ╠═902dc1c9-e1ff-4044-a01c-76141dd457d6
# ╠═346430bf-7397-426b-84a9-b322c469c89d
# ╠═31b22562-2177-4f23-9320-ce88686b0ddc
# ╠═a8e5c5df-8def-4498-95f8-2deab8b15366
# ╠═ac6dd3a5-9b1a-40d8-bb9d-0582975dc319
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

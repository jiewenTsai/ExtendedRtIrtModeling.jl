### A Pluto.jl notebook ###
# v0.19.45

using Markdown
using InteractiveUtils

# ╔═╡ 6b9551d9-840d-4afd-b304-d32413fc090d
begin
	using PlutoUI
	using LinearAlgebra
	using PolyaGammaSamplers
	using Distributions
	using MCMCChains
	using Plots
	theme(:ggplot2)
	using ProgressLogging
	using ProgressMeter
	using Random
	using Kronecker
	using StatsPlots

end

# ╔═╡ a93d8ec6-7805-11ef-002a-09f4385c2f04
md"""

# MLSA Model

In this note, I primarily collect the structs and functions used in QrSA model.

"""

# ╔═╡ 8b1cb1c8-921e-471f-a9a4-f02a161b0f11
md"""

1. 把個參數的抽樣程式打包為 `gibbsSampler()` 函數，輸入輸出學習 `jags`。應有

(Study 1)

	- [X] GibbsIrt ... Polya-Gamma Multilevel IRT Model.
	- [X] GibbsRtIrtNull

(Study 2)

	- [X] GibbsRtIrt ... 
	- [X] GibbsRtIrtQuantile ...  Bayesian Bivariate Quantile RT-IRT model.

命名原則： Gibbs + Irt/RtIrt + 擴展內容
Multilevel 是預設情況，所以 Irt 和 RtIrt 這兩個結構都沒有後綴。


用這兩個名稱建立 class，再用 sample() 函數去跑。（這個形式跟 Turing.jl 類似）。可以利用多重分派，讓 sample() 去讀不同的 struct.

這些不會用到

	- gibbsSampler2plIrt() ... 這個預設就好，不會用到。
	- gibbsSamplerMrMSa()~ ... 這個用 Lavaan 就行。
	- (null version ?) gibbsSamplerQrMSa


2. 把 input, output 整理為 struct，以便管理。應有，對應實現值（用大寫），也許用 _the_ 當前綴詞。
	- inputData ==> Data
	- inputPara ==> Para
	- outputPost ==> Post
	- outputDic ==> DIC
這幾種 struct 是通用的。所以可以先把 所有會用到的參數名稱都寫起來，再挖空即可。這些不是 funciton，是 object。

Post 的部分，`summary` 和 `mean` 另外用 `MCMCChains` 做即可。硬要放進結構中反而不順。


3. RMSE 和 Bias 只需要用到 outputMcmc 中 mean 的功能。
	- getRmse()
	- getBias()
	- getDic()


"""

# ╔═╡ 0c3642a4-c7b3-4ded-b490-77c5b73873f1
TableOfContents()

# ╔═╡ 485d8279-ee03-4cef-974e-43c145e94126
md"""
## 1.1 Struct

Struct 好像要用名詞命名，還要雙駝峰 PascalCase。因此重新想一下。

	- InputData --> Data
	- InputPara --> Para
	- OutputPost --> Post
	- OutputDic --> DIC
	- GibbsSamplerPgMlIrt --> MCMC
	- SimConditions --> Cond
	- SimSummary --> SIM

"""

# ╔═╡ ea8d8c52-e6aa-44d1-8fdb-7e8f3a4778b5
"""
    InputData --> Data

A structure to store the main components of a dataset used for modeling.

"""
struct InputData
	Y::Array
	κ::Array
	T::Array
	logT::Array
	X::Array
	function InputData(;Y=[], T=[], X=[])
		κ = Y .- 0.5
		logT = log.(T)
		return new(Y, κ, T, logT, X)
	end
end

# ╔═╡ a9602d17-d505-473e-a5f9-90e5ea792437
"""
    InputPara --> para

Represents a collection of parameter Name Lists.

This struct is designed to encapsulate various model parameters for easy management and access.

Note. Let λ = [eRa eRt]

"""
mutable struct InputPara
	ω::Array
    θ::Array
    a::Array
    b::Array
	τ::Array
	ξ::Array
	σt::Array
	λ::Array
	β::Array
	Σp::Array
	
    # 使用自動生成的構造函數，但也可以添加自定義構造函數
	function InputPara(;ω=Float64[], θ=Float64[], a=Float64[], b=Float64[], τ=Float64[], ξ=Float64[], σt=Float64[], λ=Float64[], β=Float64[], Σp=Float64[])
		return new(ω, θ, a, b, τ, ξ, σt, λ, β, Σp)
	end
end

# ╔═╡ 758449bc-5025-4193-8576-bbdeba7c66e3
"""
    OutputPostIrt --> POST
Note. For GibbsMlIrt

"""
mutable struct OutputPostMlIrt
	ra
	rt
	qr
	logLike
	mean
	function OutputPostMlIrt(Cond;
		ra=[],rt=[],qr=[],logLike=[],mean=Float64[]
	)
		ra=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain)
		qr=Array{Float64}(undef, Cond.nIter, Cond.nFeat+1, Cond.nChain)	
		logLike=Array{Float64}(undef, Cond.nIter, 1, Cond.nChain)
			
		return new(ra, rt, qr, logLike, mean)
	end
end

# ╔═╡ 3cc91952-8b1d-4e61-a120-f10c9a5989bd
"""
    OutputPost --> Post

"""
mutable struct OutputPost
	ra::Array
	rt::Array
	qr::Array
	logLike::Array
	mean
	function OutputPost(Cond; ra=[], rt=[], qr=[], logLike=[], mean=Float64[]
)
		ra=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain)
		rt=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain,)
		qr=Array{Float64}(undef, Cond.nIter, 2*(Cond.nFeat+1)+4, Cond.nChain)
		logLike=Array{Float64}(undef, Cond.nIter, 1, Cond.nChain)
		return new(ra, rt, qr, logLike, mean)
	end
end

# ╔═╡ b134cb84-4548-46c4-83e5-caa69c97697d
"""
    OutputDic --> Dic

Represents a DIC (Deviance Information Criterion) output object, used to evaluate the fit of a Bayesian model. This struct is designed to store the DIC results for model comparison and evaluation purposes.

# Fields
- `pD::Float64`: The effective number of parameters, which provides a measure of model complexity.
- `DIC::Float64`: The Deviance Information Criterion, a model selection metric that balances model fit and complexity.
"""
mutable struct OutputDic
	pD
	DIC
	function OutputDic(; pD=Float64[], DIC=Float64[])
		return new(pD, DIC)
	end
end

# ╔═╡ 7b12f2aa-a7ab-4dc3-af46-d6260a145a3b
"""
	SimConditions --> Cond
"""
struct SimConditions
	nSubj::Int
	nItem::Int
	nFeat::Int
	nIter::Int
	nChain::Int
	nBurnin::Int
	nThin::Int
	nRep::Int
	qRa::Float64
	qRt::Float64
end

# ╔═╡ 4f0acc9a-3c0f-4319-82ce-14f5e6a6a44d
md"""
### Gibbs
"""

# ╔═╡ cae950e8-75b1-4fd2-8417-7a5d5a9d054d
"""
"""
mutable struct GibbsMlIrt
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsMlIrt)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			β=randn(self.Cond.nFeat+1), 
			#Σp=[1. 0. 
			)
		return self
	end
	function GibbsMlIrt(Cond; 
	Data = [], 
	truePara = [], 
	Para = Float64[],
	Post = Float64[]	
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPostMlIrt(Cond)
		return obj
	end

end

# ╔═╡ 6c23feaa-b1a5-43f8-bd1a-423b18d9eb5f
"""
"""
mutable struct GibbsRtIrt
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsRtIrt)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			τ=randn(self.Cond.nSubj),
			ξ = zeros(self.Cond.nItem),
			σt = ones(self.Cond.nItem),
			β=randn(self.Cond.nFeat+1,2), 
			Σp=[1. 0.5; 0.5 1.] )
		return self
	end
	function GibbsRtIrt(Cond; 
	Data = [], 
	truePara = Float64[], 
	Para = Float64[],
	Post = Float64[]
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPost(Cond)
		return obj
	end

end

# ╔═╡ fc521bea-d80c-4381-9df4-3d16aeb74a11
"""
"""
mutable struct GibbsRtIrtQuantile
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsRtIrtQuantile)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			τ=randn(self.Cond.nSubj),
			ξ = zeros(self.Cond.nItem),
			σt = ones(self.Cond.nItem),
			β=randn(self.Cond.nFeat+1,2), 
			Σp=[1. 0.5; 0.5 1.] )
		return self
	end
	function GibbsRtIrtQuantile(Cond; 
	Data = Float64[], 
	truePara = Float64[], 
	Para = Float64[],
	Post = Float64[]
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPost(Cond)
		return obj
	end

end

# ╔═╡ c8685063-3ba6-42d0-af83-4b07bc573436
"""
"""
mutable struct GibbsRtIrtNull
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsRtIrtNull)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			τ=randn(self.Cond.nSubj),
			ξ=zeros(self.Cond.nItem),
			σt=ones(self.Cond.nItem),
			#β=randn(self.Cond.nFeat,2), 
			Σp=[1. 0.; 0. 1.] )
		return self
	end
	function GibbsRtIrtNull(Cond; 
	Data = Float64[], 
	truePara = Float64[], 
	Para = Float64[],
	Post = Float64[]
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPost(Cond)
		return obj
	end

end

# ╔═╡ 84dc174b-1aea-4f55-bb6f-2858a0f6d6b3
"""
"""
mutable struct SimEvaluation
	Bias
	Rmse
	Dic
	function SimEvaluation(;Bias=Float64[],Rmse=Float64[],Dic=Float64[])
		return new(Bias, Rmse, Dic)
	end
end

# ╔═╡ fac11e3a-64f7-4b0c-a3d8-a2c6bfee538b
md"""
## 1.2 Functions for preparing data
"""

# ╔═╡ 6f1f024f-0047-49b0-86af-c0ef4ceadabf
"""
"""
logistic(t) = 1 / (1 + exp(-t))

# ╔═╡ c89252d3-7f4a-471d-8d4f-2db79f9ef4f9
"""
	setCond --> Cond
"""
function setCond(; nSubj=2000, nItem=15, nFeat=3, nIter=5000, nChain=4, nBurnin=Float64[], nThin=1, nRep = 10, qRa=0.5, qRt=0.5)
	nBurnin = round(Int, nIter/2)
	return SimConditions(nSubj, nItem, nFeat, nIter, nChain, nBurnin, nThin, nRep, 	qRa, qRt)
end

# ╔═╡ 2c820be4-7dd5-44fe-8470-32a186f66bae
"""
"""
function setTrueParaRtIrt(Cond;
	a = [],
	b = [],
	β = [],
	ξ = [],
	σt = [],
	trueStdRa=1.,
	trueStdRt=sqrt(0.4),
	trueCorr=0.5
	
)
	truePara = InputPara()
	truePara.a = rand(Truncated.(Normal(1., 0.2),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., 0.5), Cond.nItem)
	truePara.β = rand(MvNormal(zeros(2), I(2)), Cond.nFeat)'
	truePara.ξ = rand(Truncated(Normal(4., 0.2 ), 0, Inf), Cond.nItem)
	truePara.σt = sqrt.(rand(Truncated.(Normal(0.3, 0.5), 0, Inf), Cond.nItem))
	trueStd = Diagonal([trueStdRa, trueStdRt])
    trueCor = [1. trueCorr; trueCorr 1.]
	truePara.Σp = trueStd * trueCor * trueStd
	
	return truePara 
end

# ╔═╡ d6cd13fe-0718-43c4-893a-98d7679497d4
"""
"""
function setTrueParaMlIrt(Cond;
	#θ = randn(COND.nSubj),
	a=[],
	b=[],
	β=[]	
)
	truePara = InputPara()
	#truePara.θ = θ
	truePara.a = rand(Truncated.(Normal(1., 0.2),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., 0.5), Cond.nItem)
	truePara.β = rand(MvNormal(zeros(1), I(1)), Cond.nFeat)'
	return truePara 
end

# ╔═╡ f2db195d-8ad1-411b-aad5-c2493fe4e431
function setData(Cond, truePara;)
	## structure
	trueX = Array{Float64}(undef, Cond.nSubj, Cond.nFeat)
	#trueX[:,1] = rand(Bernoulli(0.5), Cond.nSubj)
	trueX[:,1:end] = rand(Normal(0, 1.), Cond.nSubj, (Cond.nFeat))

	## ra
    trueMean = trueX * truePara.β
    trueStd = 1.
	truePara.θ = rand.(Normal.(trueMean, trueStd))

	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## collect
	trueData = InputData(Y=trueY, X=trueX)

	return trueData
end

# ╔═╡ 77c2580d-aa1d-436e-a9b2-41aa838ba2ac
function setDataMlIrt(Cond, truePara;)
	## structure
	trueX = Array{Float64}(undef, Cond.nSubj, Cond.nFeat)
	trueX[:,1] = rand(Bernoulli(0.5), Cond.nSubj)
	trueX[:,2:end] = rand(Normal(0, 1.), Cond.nSubj, (Cond.nFeat-1))

	## ra
    trueMean = trueX * truePara.β
    trueStd = 1.
	truePara.θ = rand.(Normal.(trueMean, trueStd))

	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## collect
	trueData = InputData(Y=trueY, X=trueX)

	return trueData
end

# ╔═╡ 5bd8ff0d-227b-489d-843c-827708126d7a
"""
"""
function setDataRtIrt(Cond, truePara)

	## structure
	trueX = Array{Float64}(undef, Cond.nSubj, Cond.nFeat)
	#trueX[:,1] = rand(Bernoulli(0.5), Cond.nSubj)
	trueX[:,1:end] = rand(Normal(0, 1.), Cond.nSubj, (Cond.nFeat))

	## ra and rt
    trueMean = trueX * truePara.β
	noise = cholesky(Symmetric(truePara.Σp)).L * randn(2,Cond.nSubj)
    trueSubj = trueMean' .+ noise
	truePara.θ = trueSubj'[:,1]
	truePara.τ = trueSubj'[:,2]
	#truePara.Σp = cov(trueSubj')
	## t(2)
	#truePara.τ = trueSubj'[:,2] ./ sqrt.(rand(Chisq(2), Cond.nSubj) ./ 2)
	## beta(0.5,0.5)
	#truePara.τ = quantile(Beta(0.5, 0.5), cdf(Normal(0, 1), trueSubj'[:,2]))
	## gamma(2.0,2.0)
	#truePara.τ = quantile(Gamma(2.0, 2.0), cdf(Normal(0, 1), trueSubj'[:,2]))


	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    μt =  truePara.ξ' .- truePara.τ 
    trueT = rand.(Truncated.(LogNormal.(μt, truePara.σt'), 0, Inf))

	## collect
	trueData = InputData(Y=trueY, X=trueX, T=trueT)


	return trueData

end

# ╔═╡ 14786bb2-28aa-498e-b7f1-eadcf4718a0e
"""
"""
function setDataRtIrtNull(Cond, truePara)

	## structure
	#trueX = Array{Float64}(undef, Cond.nSubj, Cond.nFeat)
	#trueX[:,1] = rand(Bernoulli(0.5), Cond.nSubj)
	#trueX[:,1:end] = rand(Normal(0, 1.), Cond.nSubj, (Cond.nFeat))

	## ra and rt
    #trueMean = trueX * truePara.β
	noise = cholesky(Symmetric(truePara.Σp)).L * randn(2,Cond.nSubj)
    trueSubj = noise
	truePara.θ = trueSubj'[:,1]
	truePara.τ = trueSubj'[:,2]
	#truePara.Σp = cov(trueSubj')
	## t(2)
	#truePara.τ = trueSubj'[:,2] ./ sqrt.(rand(Chisq(2), Cond.nSubj) ./ 2)
	## beta(0.5,0.5)
	#truePara.τ = quantile(Beta(0.5, 0.5), cdf(Normal(0, 1), trueSubj'[:,2]))
	## gamma(2.0,2.0)
	#truePara.τ = quantile(Gamma(2.0, 2.0), cdf(Normal(0, 1), trueSubj'[:,2]))


	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    μt =  truePara.ξ' .- truePara.τ 
    trueT = rand.(Truncated.(LogNormal.(μt, truePara.σt'), 0, Inf))

	## collect
	trueData = InputData(Y=trueY, T=trueT)


	return trueData

end

# ╔═╡ cdc7fbd2-8315-4c50-832f-0e1e353f1172
md"""
## 1.2 _draw_-Functions
"""

# ╔═╡ 5a502f85-2624-4ba9-abea-7428f1989816
md" ### 1.2.1 2PL IRT"

# ╔═╡ 1000ab97-3d0b-4f3f-984a-c7342d0e79d8
"""
	drawRaPgRandomVariable(Para) --> ω
"""
function drawRaPgRandomVariable(Para)
    #η = Para.θ * Para.a' .- Para.b'
	η = Para.a' .* (Para.θ .- Para.b')
	ω = rand.(PolyaGammaPSWSampler.(1, η))
    return ω
end

# ╔═╡ fbed8560-7adc-4ad4-a14f-4dbda25ad67e
"""
"""
function drawSubjAbilityNull(Cond,Data,Para)
    #x = [ones(Cond.nSubj) Data.X]

	θμ₀ = 0.
	θσ₀² = 1. #Para.Σp[1,1]

    parV = 1 ./(1 ./θσ₀² .+ sum(Para.a'.^2 .* Para.ω, dims=2))
	parM = parV .* (θμ₀ ./θσ₀² .+ sum( Para.a' .* (Data.κ .+ Para.a' .* Para.b' .* Para.ω ), dims=2 ))
	#θ = rand.(Normal.(parM, sqrt.(parV))) 
    θ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), -10,10))
	return θ
end

# ╔═╡ 970a24bf-8ce5-4e9e-89a9-fc7a6c24d36b
"""
"""
function drawSubjAbility(Cond,Data,Para)
    x = [ones(Cond.nSubj) Data.X]

	θμ₀ = x * Para.β[:,1] 
	θσ₀² = 1. #Para.Σp[1,1]

    parV = 1 ./(1 ./θσ₀² .+ sum(Para.a'.^2 .* Para.ω, dims=2))
	parM = parV .* (θμ₀ ./θσ₀² .+ sum( Para.a' .* (Data.κ .+ Para.a' .* Para.b' .* Para.ω ), dims=2 ))
	#θ = rand.(Normal.(parM, sqrt.(parV))) 
    θ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), -10,10))
	return θ
end

# ╔═╡ 0f6abf0d-aab8-4b40-806f-978f9cbb2f37
"""
"""
function drawItemDiscrimination(Data, Para; μa₀=1., σa₀=1.)
	parV = 1 ./ (1/σa₀^2 .+ sum( (Para.θ .- Para.b').^2 .*  Para.ω, dims=1))
    parM = parV .* ( μa₀/σa₀^2  .+  sum( Data.κ .* (Para.θ .- Para.b'), dims=1) )
    a = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf)  ) 
    return a'
end

# ╔═╡ 9c437f8f-d1bc-4218-9f83-5a438c8c22b6
"""
"""
function drawItemDifficulty(Data, Para; μb₀=0., σb₀=1e+10)
    parV = 1 ./ (1/σb₀^2 .+  sum( Para.a'.^2 .* Para.ω, dims=1)')
    parM = parV .* ( μb₀/σb₀^2 .- Para.a' .* sum( Data.κ .- Para.θ*Para.a'.*Para.ω, dims=1) )'
    b = rand.(Normal.(parM, sqrt.(parV))) 
    return b
end

# ╔═╡ a975df83-7e0f-45eb-bb1b-c2c673d5d6e9
md"""
### 1.2.2 RT Model
"""

# ╔═╡ de279a3f-0c73-4542-8aaa-9d397d2c8196
"""
"""
function drawSubjSpeedNull(Cond,Data,Para )
	#x = [ones(Cond.nSubj) Data.X]

	τμ₀ = 0. #x * Para.β[:,2] 
	τσ₀² = 1. #Para.Σp[2,2] 
	
    parV = 1 ./ (1 ./ τσ₀² .+ sum( 1. ./  Para.σt.^2, dims=1) )
    parM = parV .* (τμ₀ ./ τσ₀² .+ sum( (Para.ξ' .- Data.logT) ./ Para.σt'.^2, dims=2))  
    #τ = rand.(Normal.(parM, sqrt.(parV)))
    τ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), -10,10))
    return τ
end

# ╔═╡ 8ba1f85f-7498-4b3e-82e2-36bf29998265
"""
"""
function drawSubjSpeed(Cond,Data,Para )
	x = [ones(Cond.nSubj) Data.X]

	τμ₀ = x * Para.β[:,2] 
	τσ₀² = Para.Σp[2,2] 
	
    parV = 1 ./ (1 ./ τσ₀² .+ sum( 1. ./  Para.σt.^2, dims=1) )
    parM = parV .* (τμ₀ ./ τσ₀² .+ sum( (Para.ξ' .- Data.logT) ./ Para.σt'.^2, dims=2))  
    #τ = rand.(Normal.(parM, sqrt.(parV)))
    τ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), -10,10))
    return τ
end

# ╔═╡ 718c896b-5746-4ad9-9e0b-fa7af2dc56cd
"""
"""
function drawItemIntensity(Cond,Data,Para; μξ=mean(Data.logT), σξ=1e+10)
    parV = 1 ./(1/σξ^2 .+ Cond.nSubj ./ Para.σt.^2)
    parM = parV .* (μξ/σξ^2 .+ sum(Data.logT .+ Para.τ, dims=1) ./ Para.σt'.^2)'
    ξ = rand.(Truncated.(Normal.(parM, sqrt.(parV)), 0, Inf))
    return ξ
end

# ╔═╡ c7c96185-af24-4c38-b7af-c1897acd99bd
function drawItemTimeResidual(Cond,Data,Para;δa=1e-10, δb=1e-10)
    parA = δa + Cond.nSubj/2
    parB = δb .+ sum( (Data.logT .- Para.ξ' .+ Para.τ).^2, dims=1)./2
    σt = rand.(InverseGamma.(parA, parB'))
    return sqrt.(σt)
end

# ╔═╡ 888f5005-a524-42e6-b8a2-879418e4248a
md"""
### 1.2.3 Structure
"""

# ╔═╡ 89f381ed-f0fb-462e-b014-805c1575e139

"""
    drawQrWeights(θ, τ, x; qRa=0.5, qRt=0.5, β, Σp) --> λ
    (New!)
"""
function drawQrWeights(Cond,Data,Para)
    x = [ones(Cond.nSubj) Data.X]

    k1Ra = (1 - 2 * Cond.qRa) / (Cond.qRa * (1 - Cond.qRa))
    k2Ra = 2 / (Cond.qRa * (1 - Cond.qRa))
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
	
    parRaM = sqrt.(2 * k2Ra .+ k1Ra.^2) ./ abs.(Para.θ .- x* Para.β[:,1])
    parRaL = (2 * k2Ra .+ k1Ra.^2) ./ (k2Ra * 1 ./ Para.Σp[1,1])

    parRtM = sqrt.(2 * k2Rt .+ k1Rt.^2) ./ abs.(Para.τ .- x * Para.β[:,2])
    parRtL = (2 * k2Rt .+ k1Rt.^2) ./ (k2Rt * 1 ./ Para.Σp[2,2])

    eRa =  1 ./ rand.(InverseGaussian.(parRaM, parRaL))
    eRt =  1 ./ rand.(InverseGaussian.(parRtM, parRtL))


    return [eRa eRt]
end

# ╔═╡ 1f9287ca-cb85-4e7a-8e5f-32a343f94e2c
"""
"""
function getSubjCoefficientsMlIrt(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )
	η = Para.θ
	x = [ones(Cond.nSubj) Data.X]

    #invΩ = inv(Para.Σp) 
    #parV = inv( 1/σβ₀^2 .+ invΩ .* x'x)  # Regularization with identity matrix
    #parM = parV * (μβ₀/σβ₀^2 .+ vec(x'* Para.θ * invΩ') )

	β = x'x \ x'*η

    ## Generate random coefficients
    #β = parM .+ cholesky(Symmetric(parV)).L * randn((Cond.nFeat))

    return β
end

# ╔═╡ 9a6816dd-fcd4-4fe8-a8e0-f17038976077
"""
"""
function drawSubjCoefficients(Cond,Data,Para; μβ₀=0., σβ₀=1e+10 )

    η = [Para.θ Para.τ]
    x = [ones(Cond.nSubj) Data.X]

    invΩ = inv(Symmetric(Para.Σp))
	parV = inv( 1/σβ₀^2 .+  invΩ ⊗ x'x) 
    parM = parV * ( μβ₀/σβ₀^2 .+ vec(x'* η * invΩ') )

    ## Generate random coefficients
    β = parM .+ cholesky(Symmetric(parV)).L * randn((Cond.nFeat+1)*2)

    return reshape(β, Cond.nFeat+1,2)
end

# ╔═╡ 94583a08-4dfb-4579-b8a2-539833c6645d
"""
"""
function drawSubjCoefficientsQr(Cond, Data, Para; βμ=0., βσ=1e+10 )

    η = [Para.θ Para.τ]
    x = [ones(Cond.nSubj) Data.X]
	eRa = Para.λ[:,1]
	eRt = Para.λ[:,2]
    k1Ra = (1 - 2 * Cond.qRa) / (Cond.qRa * (1 - Cond.qRa))
    k2Ra = 2 / (Cond.qRa * (1 - Cond.qRa))
    k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
    k1e = [k1Ra.* eRa k1Rt.* eRt]

    invΩRa = sqrt.(sum(k2Ra*eRa,dims=1))
    invΩRt = sqrt.(sum(k2Rt*eRt,dims=1))
    sqrtK2e = [invΩRa 0.; 0. invΩRt]


    #invΩ = inv(Para.Σp) 
    invΩ = sqrtK2e * inv(Para.Σp) * sqrtK2e
	parV = inv(invΩ ⊗ x'x) 
    parM = parV * (vec(x'*(η .- k1e) * invΩ') )
	#parM = parV * (vec(x'*(η .- k1e)) )

	#parV = inv( 1/βσ^2 .+  invΩ ⊗ x'x)  # Regularization with identity matrix
    #parM = parV * ( βμ/βσ^2 .+ vec(x'*(η .- k1e) * invΩ') )

    β = parM .+ cholesky(Symmetric(parV)).L * randn((Cond.nFeat+1)*2)

    return reshape(β, Cond.nFeat+1,2)
end

# ╔═╡ 26771447-5b21-4733-95d5-38609d0239d8
"""
"""
function getSubjCoefficients(Cond, Data, Para )
	η = [Para.θ Para.τ]
	x = [ones(Cond.nSubj) Data.X]
	invΩ = inv(Para.Σp)
	#β = x'x \ x'*η
	β = (invΩ ⊗ x'x) \ vec(x'*η * invΩ')
	#β = (invΩ ⊗ x'x) \ vec(x'*η )
	#β = Data.X'Data.X \ (Data.X' * η * inv(Para.Σp))

	β = reshape(β, Cond.nFeat+1,2)
	return β
end

# ╔═╡ 07522aa2-ce40-4a90-877c-ff1b5fadd789
"""
    drawSubjCovariance(;θ, τ, nSubj)
    (New)
"""
function drawSubjCovariance(Cond, Data, Para)
    η = [Para.θ Para.τ]
	#meanη = [mean(Para.θ) mean(Para.τ)]
    x = [ones(Cond.nSubj) Data.X]
	e = η .- x * Para.β 

    ee = e'e
    s =  rand(InverseWishart(Cond.nSubj+3, ee + I(2) ))
	
	#s[1,2] = s[2,1] = s[1,2] ./ sqrt(xi[1] * xi[2])
    #s[1,2] = s[2,1] = s[1,2] ./ sqrt.(s[1,1] * s[2,2])
	#s[1,1] = s[1,1] ./ s[1,1]


	L = cholesky(s).L
	L[2,1] = L[2,1] ./ (L[1,1] * L[2,2])

	ss = L*L'

    return ss
end

# ╔═╡ c5c72069-ab2b-4f79-a80b-d256441a8a6c
"""
    drawSubjCovariance(;θ, τ, nSubj)
    (New)
"""
function drawSubjCovarianceNull(Cond, Data, Para)
    η = [Para.θ Para.τ]
	#meanη = [mean(Para.θ) mean(Para.τ)]
	#e = η .- meanη
    ee = η'η
	#s = e'e ./ Cond.nSubj
	s =  rand(InverseWishart(Cond.nSubj+3, ee + I(2)))

	#s[1,2] = s[2,1] = s[1,2] ./ sqrt.(s[1,1] * s[2,2])
	#s[1,1] = s[1,1] ./ s[1,1]
	#s[1,2] =  s[1,2] ./ sqrt.(s[1,1] * s[2,2])
	#s[2,1] =  s[2,1] ./ sqrt.(s[1,1] * s[2,2])
    #ss = diagm([s[1,1].^-0.5,1.]) * s * diagm([s[1,1].^-0.5, 1.]) 

    return s #Σp
end

# ╔═╡ 6c19bce2-70f2-43d9-9fdf-8cf1bb505ac1
md"""

## ⭐ 打包區

"""

# ╔═╡ 30f39851-d610-4faa-adf8-1b81a42c09a5
md"""
## 1.3 *get-*Functions.

這邊凡事要「取得一個值」的函數，都以 *get* 當函數名稱。

"""

# ╔═╡ e1c7895c-a770-4e97-bb21-04a34baaa14c
"""
"""
getRmse(a,b) = mean(sqrt.(mean.((a - b).^2)))

# ╔═╡ 4626846f-4c5a-4171-bcec-9d64bcc07bac
"""
"""
getBias(a,b) = mean(mean.(a - b))

# ╔═╡ c7edc678-e100-4bfb-a8e8-55486ba8d225
function evaluate(MCMC::GibbsRtIrt)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b) getRmse(MCMC.Post.mean.β, MCMC.truePara.β)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b) getBias(MCMC.Post.mean.β, MCMC.truePara.β)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end

# ╔═╡ 37e82d85-bfeb-482c-a8e4-7ccd1e34431d
function evaluate(MCMC::GibbsRtIrtQuantile)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b) getRmse(MCMC.Post.mean.β, MCMC.truePara.β)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b) getBias(MCMC.Post.mean.β, MCMC.truePara.β)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end

# ╔═╡ 929f5926-a6be-4db9-b5d7-69154409537f
function evaluate(MCMC::GibbsRtIrtNull)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end

# ╔═╡ aaa4e2de-90c9-43c3-b36d-8eb23842419f
"""
"""
function getLogLikelihood2plIrt(Data,Post)
    pr = logistic.(Post.mean.a' .* (Post.mean.θ .- Post.mean.b'))
    sum(log.(pr).*Data.Y + log.(1 .- pr).*(1 .-Data.Y))
end

# ╔═╡ 8a383a20-8a9b-428b-8e11-d8e3229bd3ff
"""
	getLogLikelihood(MCMC,Data)
"""
function getLogLikelihoodMlIrt(Cond,Data; P=Para)
	x = [ones(Cond.nSubj) Data.X]
    pr = P.a' .* (P.θ .- P.b')
    #μ = ξ' .- τ
    #μOfθ = [ones(Cond.nSubj) Data.X] * P.β
	μOfθ = x * P.β
    logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μOfθ, 1.), P.θ))]

    return sum(logPdf)
end

# ╔═╡ 91235ba9-0fb9-4452-86a3-d4a88db41b7e
"""
	sample!(MCMC::GibbsMlIrt)
"""
function sample!(MCMC::GibbsMlIrt; intercept=false, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end

	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain
	
		## structural
		Para.β = getSubjCoefficientsMlIrt(Cond,Data,Para)
		if intercept==false
			Para.β[1] = 0.
		end
		##Para.Σp = 1. #drawSubjVariance(COND,DATA,PARA;)
		
	
		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.a = drawItemDiscrimination(Data, Para)
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		Para.b = drawItemDifficulty(Data,Para)
		Para.θ = drawSubjAbility(Cond, Data, Para)

		## collect
		Post.qr[m,:,l] = vec(Para.β)
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]

		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodMlIrt(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(
		θ = mean(Post.ra[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		a = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		b = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec,
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:(1*(Cond.nFeat+1)), :], dims=(1,3)) |> vec
	)

	return MCMC
end

# ╔═╡ 762a1217-6892-4b1e-9ee4-ed6de7f63f7e
"""
"""
function getLogLikelihoodRtIrt(Cond,Data; P=Para)
	x = [ones(Cond.nSubj) Data.X]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.ξ' .- P.τ
	η = [P.θ P.τ]
    μη = x * reshape(P.β, Cond.nFeat+1, 2)
	Ση = reshape(P.Σp, 2,2) 
	logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, P.σt'), Data.logT)) sum([logpdf(MvNormal(μη[i,:], Ση), η[i,:]) for i in 1:Cond.nSubj])]

    return sum(logPdf)
end


# ╔═╡ ad33ce91-0f89-4520-89fb-3f06645d4df3
"""
	sample!(MCMC::GibbsRtIrtQuantile)
"""
function sample!(MCMC::GibbsRtIrtQuantile; intercept=false, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	
		## structural
		Para.λ = drawQrWeights(Cond,Data,Para)
		Para.β = drawSubjCoefficientsQr(Cond,Data,Para)
		if intercept==false
			Para.β[1,:] = [0. 0.]
		end
		Para.Σp = drawSubjCovariance(Cond, Data, Para)
		
		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.b = drawItemDifficulty(Data,Para)
		Para.a = drawItemDiscrimination(Data, Para)
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		Para.θ = drawSubjAbility(Cond, Data, Para)

		
		
		## response time part
		Para.ξ = drawItemIntensity(Cond,Data,Para)
	    Para.σt = drawItemTimeResidual(Cond,Data,Para)
		Para.τ = drawSubjSpeed(Cond,Data,Para)

		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.τ; Para.ξ; Para.σt]
		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodRtIrt(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(

		## structure
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:(2*(Cond.nFeat+1)), :], dims=(1,3)) |> vec,
		Σp = mean(Post.qr[(Cond.nBurnin+1):end, (2*(Cond.nFeat+1)+1):end, :], dims=(1,3)) |> vec,

		## ra
		θ = mean(Post.ra[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		a = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		b = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec,

		## rt
		τ = mean(Post.rt[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		ξ = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		σt = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec

	)

	return MCMC
end

# ╔═╡ 2c4f365d-321b-4588-a49d-5a29dc4432b0
"""
	sample(MCMC::GibbsRtIrt)
"""
function sample!(MCMC::GibbsRtIrt; intercept=false, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	

		## structural
		#Para.λ = drawQrWeights(Cond,Data,Para)
		Para.β = drawSubjCoefficients(Cond,Data,Para)
		if intercept==false
			Para.β[1,:] = [0. 0.]
		end
		Para.Σp = drawSubjCovariance(Cond, Data, Para)
		
		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.b = drawItemDifficulty(Data,Para)
		Para.a = drawItemDiscrimination(Data, Para)
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		Para.θ = drawSubjAbility(Cond, Data, Para)
		
		
		## response time part
	    Para.ξ = drawItemIntensity(Cond,Data,Para)
	    Para.σt = drawItemTimeResidual(Cond,Data,Para)
		Para.τ = drawSubjSpeed(Cond,Data,Para)



		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.τ; Para.ξ; Para.σt]
		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodRtIrt(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(

		## structure
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:(2*(Cond.nFeat+1)), :], dims=(1,3)) |> vec,
		Σp = mean(Post.qr[(Cond.nBurnin+1):end, (2*(Cond.nFeat+1)+1):end, :], dims=(1,3)) |> vec,

		## ra
		θ = mean(Post.ra[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		a = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		b = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec,

		## rt
		τ = mean(Post.rt[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		ξ = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		σt = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec

	)

	return MCMC
end

# ╔═╡ 335b4fc4-a079-4107-8726-2fa802468041
"""
"""
function getLogLikelihoodRtIrtNull(Cond,Data; P=Para)
	#x = [ones(Cond.nSubj) Data.X]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.ξ' .- P.τ
	η = [P.θ P.τ]
    μη = zeros(Cond.nSubj, 2) #x * reshape(P.β, Cond.nFeat+1, 2)
	Ση = reshape(P.Σp, 2,2)
	logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, P.σt'), Data.logT)) sum([logpdf(MvNormal(μη[i,:], Ση), η[i,:]) for i in 1:Cond.nSubj])]

    return sum(logPdf)
end


# ╔═╡ ea14de6e-15b0-4c3f-af7b-cdf57de0bf2e
"""
	sample(MCMC::GibbsRtIrtNull)
"""
function sample!(MCMC::GibbsRtIrtNull, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	
		## structural
		#Para.λ = drawQrWeights(Cond,Data,Para)
		Para.β = zeros(Cond.nFeat+1, 2)
		Para.Σp = drawSubjCovarianceNull(Cond, Data, Para)
		
		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.b = drawItemDifficulty(Data,Para)
		Para.a = drawItemDiscrimination(Data, Para)	
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		
		Para.θ = drawSubjAbilityNull(Cond, Data, Para)
		
		## response time part
	    Para.ξ = drawItemIntensity(Cond,Data,Para)
	    Para.σt = drawItemTimeResidual(Cond,Data,Para)
		Para.τ = drawSubjSpeedNull(Cond,Data,Para)

		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.τ; Para.ξ; Para.σt]
		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodRtIrtNull(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(

		## structure
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:(2*(Cond.nFeat+1)), :], dims=(1,3)) |> vec,
		Σp = mean(Post.qr[(Cond.nBurnin+1):end, (2*(Cond.nFeat+1)+1):end, :], dims=(1,3)) |> vec,

		## ra
		θ = mean(Post.ra[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		a = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		b = mean(Post.ra[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec,

		## rt
		τ = mean(Post.rt[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		ξ = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		σt = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec

	)

	return MCMC
end

# ╔═╡ 5eeb8329-0c72-49b7-9179-07df457861fc
md"""
	- println('DIC info (using the rule, pD = D̄ - D̂)')
	- println(pD = $(round(DIC.pD, digits=3)) and DIC = $(round(DIC.DIC, digits=3)))
	- println('DIC is an estimate of expected predictive error (lower deviance is better).')

"""

# ╔═╡ e5c47428-e543-4836-8f93-57ea61aec05d
"""
	getDicMlIrt(Cond, Data, Post)
"""
function getDic(MCMC::GibbsMlIrt)
	DIC = OutputDic()
	
	## compute dic
	D̂ = -2*getLogLikelihoodMlIrt(MCMC.Cond, MCMC.Data,P=MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD


	return DIC
	
end

# ╔═╡ 39ea5fbf-3d06-4050-b8b1-12c9eddff127
"""
"""
function getDic(MCMC::GibbsRtIrt)
	DIC = OutputDic()	
	## compute dic
	D̂ = -2*getLogLikelihoodRtIrt(MCMC.Cond, MCMC.Data, P= MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD
	return DIC
	
end

# ╔═╡ 0c3be9d9-89d1-4184-b80a-adbe8c86d15c
"""
"""
function getDic(MCMC::GibbsRtIrtQuantile)
	DIC = OutputDic()	
	## compute dic
	D̂ = -2*getLogLikelihoodRtIrt(MCMC.Cond, MCMC.Data, P= MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD
	return DIC
	
end

# ╔═╡ 563c6957-d66c-4b74-8442-ad0c74aeff9f
"""
"""
function getDic(MCMC::GibbsRtIrtNull)
	DIC = OutputDic()	
	## compute dic
	D̂ = -2*getLogLikelihoodRtIrtNull(MCMC.Cond, MCMC.Data, P= MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD
	return DIC
	
end

# ╔═╡ 3fe057b5-7dda-4843-8f58-0cde265cfc80
md"### Testing "

# ╔═╡ d4b00edb-7a00-480f-a8dd-280daf5087e1
"""
begin
	Cond2 = setCond(qRa=0.05, qRt=0.95, nChain=3, nIter=3000)
	truePara2 = setTrueParaSa(Cond2, trueCorr=0.8)
	Data2 = setDataRtIrt(Cond2, truePara2)
	Cond3 = setCond(qRa=0.5, qRt=0.5, nChain=3, nIter=3000)
end
"""

# ╔═╡ 6a8d0b9b-c11d-4f79-b2aa-80d47b3b8094
#using CSV, DataFrames

# ╔═╡ 1af30d48-3553-4fe1-8a2a-b5c598740b71
#Data0929 = CSV.read("Data0929.csv", DataFrame)

# ╔═╡ 039c4bba-e87f-430c-bbe1-8e6c2d98df02
"""
Data3 = InputData(
    Y=Matrix(Data0929[:,1:15]),
    T=exp.(Matrix(Data0929[:,16:30])),
    X=Matrix(Data0929[:,31:33])
)
"""

# ╔═╡ 38b93062-57e4-486b-af9c-aa4a73080f01
#MCMC = GibbsRtIrt(Cond3, Data=Data2, truePara=truePara2)

# ╔═╡ 49b24a83-2385-4f9d-a607-fd685a953d89
#sample!(MCMC)

# ╔═╡ b9eb3646-1f58-4a14-8615-c5850a31f65d
#MCMC.Post.mean

# ╔═╡ 2ba37c20-f1e1-4875-b9ec-6ba0647f5b1b
#MCMC.truePara.Σp

# ╔═╡ a59e9814-60e9-4c91-ac22-3883eb154c9d
#getRmse(MCMC.truePara.a, MCMC.Post.mean.a)

# ╔═╡ 8259b697-54e3-448a-9332-4ad9e7c7fc45
#getRmse(MCMC.truePara.θ, MCMC.Post.mean.θ)

# ╔═╡ 5d1d2385-bba0-4ac4-8336-938078c1a7e8
#getRmse(MCMC.truePara.τ, MCMC.Post.mean.τ)

# ╔═╡ b98d9dce-6c4f-4a6f-b4d4-19b42bf96c08
#marginalhist(MCMC.Post.mean.θ, MCMC.Post.mean.τ)

# ╔═╡ 14a33fb9-0733-4db8-90b9-89c370458b98
#Cond8 = setCond(nSubj=1000, nItem=15, qRa=0.5, qRt=0.5, nChain=3, nIter=3000)

# ╔═╡ 5cdb09bd-8e35-4b46-a14f-b32de5ac0429
#truePara8 = setTrueParaMlIrt(Cond8)

# ╔═╡ 76767de6-85d2-4af7-9c54-3c4f0070002d
#Data8 = setDataMlIrt(Cond8, truePara8)

# ╔═╡ e4564810-15fc-4820-9da1-b0818bd6be0e
#MCMC8 = GibbsMlIrt(Cond8, Data=Data8, truePara=truePara8)

# ╔═╡ 4c88ac80-d6a4-4192-a4e8-e4bc23fecc59
#sample!(MCMC8, intercept=true, itemtype="1pl")

# ╔═╡ 708eb634-9e44-4164-9afd-33225c8a9013
#MCMC8.Post.mean

# ╔═╡ fcba6869-3224-46e3-9b3b-b754e7671499
#plot(MCMC8.Post.mean.θ)

# ╔═╡ 44c6092d-9a02-4718-8e44-db407b5f124e
#MCMC8.Post.logLike

# ╔═╡ 8c3cee3e-e283-49d2-9425-321c128d528a
#getDic(MCMC8)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Kronecker = "2c470bb0-bcc8-11e8-3dad-c9649493f05e"
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
MCMCChains = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PolyaGammaSamplers = "99ff7fc7-8a0f-4729-8284-81f1989d3fc6"
ProgressLogging = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"
Random = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
StatsPlots = "f3b207a7-027a-5e70-b257-86293d7955fd"

[compat]
Distributions = "~0.25.109"
Kronecker = "~0.5.5"
MCMCChains = "~5.7.1"
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
PolyaGammaSamplers = "~0.1.0"
ProgressLogging = "~0.1.4"
ProgressMeter = "~1.10.2"
StatsPlots = "~0.15.7"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.5"
manifest_format = "2.0"
project_hash = "fde1fe7d62d1c03c339d052cbada2db27120d3df"

[[deps.AbstractFFTs]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d92ad398961a3ed262d8bf04a1a2b8340f915fef"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.5.0"
weakdeps = ["ChainRulesCore", "Test"]

    [deps.AbstractFFTs.extensions]
    AbstractFFTsChainRulesCoreExt = "ChainRulesCore"
    AbstractFFTsTestExt = "Test"

[[deps.AbstractMCMC]]
deps = ["BangBang", "ConsoleProgressMonitor", "Distributed", "LogDensityProblems", "Logging", "LoggingExtras", "ProgressLogging", "Random", "StatsBase", "TerminalLoggers", "Transducers"]
git-tree-sha1 = "87e63dcb990029346b091b170252f3c416568afc"
uuid = "80f14c24-f653-4e6a-9b94-39d6b0f70001"
version = "4.4.2"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.AbstractTrees]]
git-tree-sha1 = "2d9c9a55f9c93e8887ad391fbae72f8ef55e1177"
uuid = "1520ce14-60c1-5f80-bbc7-55ef81b5835c"
version = "0.4.5"

[[deps.Adapt]]
deps = ["LinearAlgebra", "Requires"]
git-tree-sha1 = "6a55b747d1812e699320963ffde36f1ebdda4099"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "4.0.4"
weakdeps = ["StaticArrays"]

    [deps.Adapt.extensions]
    AdaptStaticArraysExt = "StaticArrays"

[[deps.AliasTables]]
deps = ["PtrArrays", "Random"]
git-tree-sha1 = "9876e1e164b144ca45e9e3198d0b689cadfed9ff"
uuid = "66dad0bd-aa9a-41b7-9441-69ab47430ed8"
version = "1.1.3"

[[deps.ArgCheck]]
git-tree-sha1 = "a3a402a35a2f7e0b87828ccabbd5ebfbebe356b4"
uuid = "dce04be8-c92d-5529-be00-80e4d2c0e197"
version = "2.3.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Arpack]]
deps = ["Arpack_jll", "Libdl", "LinearAlgebra", "Logging"]
git-tree-sha1 = "9b9b347613394885fd1c8c7729bfc60528faa436"
uuid = "7d9fca2a-8960-54d3-9f78-7d1dccf2cb97"
version = "0.5.4"

[[deps.Arpack_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "OpenBLAS_jll", "Pkg"]
git-tree-sha1 = "5ba6c757e8feccf03a1554dfaf3e26b3cfc7fd5e"
uuid = "68821587-b530-5797-8361-c406ea357684"
version = "3.5.1+1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "01b8ccb13d68535d73d2b0c23e39bd23155fb712"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.1.0"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "16351be62963a67ac4083f748fdb3cca58bfd52f"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.7"

[[deps.BangBang]]
deps = ["Compat", "ConstructionBase", "InitialValues", "LinearAlgebra", "Requires", "Setfield", "Tables"]
git-tree-sha1 = "7aa7ad1682f3d5754e3491bb59b8103cae28e3a3"
uuid = "198e06fe-97b7-11e9-32a5-e1d131e6ad66"
version = "0.3.40"

    [deps.BangBang.extensions]
    BangBangChainRulesCoreExt = "ChainRulesCore"
    BangBangDataFramesExt = "DataFrames"
    BangBangStaticArraysExt = "StaticArrays"
    BangBangStructArraysExt = "StructArrays"
    BangBangTypedTablesExt = "TypedTables"

    [deps.BangBang.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    StaticArrays = "90137ffa-7385-5640-81b9-e52037218182"
    StructArrays = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
    TypedTables = "9d95f2ec-7b3d-5a63-8d20-e2491e220bb9"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Baselet]]
git-tree-sha1 = "aebf55e6d7795e02ca500a689d326ac979aaf89e"
uuid = "9718e550-a3fa-408a-8086-8db961cd8217"
version = "0.1.1"

[[deps.BitFlags]]
git-tree-sha1 = "0691e34b3bb8be9307330f88d1a3c3f25466c24d"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.9"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a2f1c8c668c8e3cb4cca4e57a8efdb09067bb3fd"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+2"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra"]
git-tree-sha1 = "71acdbf594aab5bbb2cec89b208c41b4c411e49f"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.24.0"
weakdeps = ["SparseArrays"]

    [deps.ChainRulesCore.extensions]
    ChainRulesCoreSparseArraysExt = "SparseArrays"

[[deps.Clustering]]
deps = ["Distances", "LinearAlgebra", "NearestNeighbors", "Printf", "Random", "SparseArrays", "Statistics", "StatsBase"]
git-tree-sha1 = "9ebb045901e9bbf58767a9f34ff89831ed711aae"
uuid = "aaaa29a8-35af-508c-8bc3-b662a17a0fe5"
version = "0.15.7"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "b8fe8546d52ca154ac556809e10c75e6e7430ac8"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.5"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "b5278586822443594ff615963b0c09755771b3e0"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.26.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"
weakdeps = ["SpecialFunctions"]

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "362a287c3aa50601b0bc359053d5c2468f0e7ce0"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.11"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "8ae8d32e09f0dcf42a36b90d4e17f5dd2e4c4215"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.16.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.1.1+0"

[[deps.CompositionsBase]]
git-tree-sha1 = "802bb88cd69dfd1509f6670416bd4434015693ad"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.2"
weakdeps = ["InverseFunctions"]

    [deps.CompositionsBase.extensions]
    CompositionsBaseInverseFunctionsExt = "InverseFunctions"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "ea32b83ca4fefa1768dc84e504cc0a94fb1ab8d1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.2"

[[deps.ConsoleProgressMonitor]]
deps = ["Logging", "ProgressMeter"]
git-tree-sha1 = "3ab7b2136722890b9af903859afcf457fa3059e8"
uuid = "88cd18e8-d9cc-4ea6-8889-5259c0d15c8b"
version = "0.1.2"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "d8a9c0b6ac2d9081bf76324b39c78ca3ce4f0c98"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.5.6"
weakdeps = ["IntervalSets", "StaticArrays"]

    [deps.ConstructionBase.extensions]
    ConstructionBaseIntervalSetsExt = "IntervalSets"
    ConstructionBaseStaticArraysExt = "StaticArrays"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Dbus_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fc173b380865f70627d7dd1190dc2fce6cc105af"
uuid = "ee1fde0b-3d02-5ea6-8484-8dfef6360eab"
version = "1.14.10+0"

[[deps.DefineSingletons]]
git-tree-sha1 = "0fba8b706d0178b4dc7fd44a96a92382c9065c2c"
uuid = "244e2a9f-e319-4986-a169-4d1fe445cd52"
version = "0.1.2"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.Distances]]
deps = ["LinearAlgebra", "Statistics", "StatsAPI"]
git-tree-sha1 = "66c4c81f259586e8f002eacebc177e1fb06363b0"
uuid = "b4f34e82-e78d-54a5-968a-f98e89d6e8f7"
version = "0.10.11"
weakdeps = ["ChainRulesCore", "SparseArrays"]

    [deps.Distances.extensions]
    DistancesChainRulesCoreExt = "ChainRulesCore"
    DistancesSparseArraysExt = "SparseArrays"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.Distributions]]
deps = ["AliasTables", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SpecialFunctions", "Statistics", "StatsAPI", "StatsBase", "StatsFuns"]
git-tree-sha1 = "9c405847cc7ecda2dc921ccf18b47ca150d7317e"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.109"

    [deps.Distributions.extensions]
    DistributionsChainRulesCoreExt = "ChainRulesCore"
    DistributionsDensityInterfaceExt = "DensityInterface"
    DistributionsTestExt = "Test"

    [deps.Distributions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    DensityInterface = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
    Test = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1c6317308b9dc757616f0b5cb379db10494443a7"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.6.2+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FFTW]]
deps = ["AbstractFFTs", "FFTW_jll", "LinearAlgebra", "MKL_jll", "Preferences", "Reexport"]
git-tree-sha1 = "4820348781ae578893311153d69049a93d05f39d"
uuid = "7a1cc6ca-52ef-59f5-83cd-3a7055c09341"
version = "1.8.0"

[[deps.FFTW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c6033cc3892d0ef5bb9cd29b7f2f0331ea5184ea"
uuid = "f5851436-0d7a-5f13-b9de-f02708fd171a"
version = "3.3.10+0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FillArrays]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "0653c0a2396a6da5bc4766c43041ef5fd3efbe57"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "1.11.0"
weakdeps = ["PDMats", "SparseArrays", "Statistics"]

    [deps.FillArrays.extensions]
    FillArraysPDMatsExt = "PDMats"
    FillArraysSparseArraysExt = "SparseArrays"
    FillArraysStatisticsExt = "Statistics"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "05882d6995ae5c12bb5f36dd2ed3f61c98cbb172"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.5"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Zlib_jll"]
git-tree-sha1 = "db16beca600632c95fc8aca29890d83788dd8b23"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.96+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.Formatting]]
deps = ["Logging", "Printf"]
git-tree-sha1 = "fb409abab2caf118986fc597ba84b50cbaf00b87"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.3"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "5c1d8ae0efc6c2e7b1fc502cbe25def8f661b7bc"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.2+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1ed150b39aebcc805c26b93a8d0122c940f64ce2"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.14+0"

[[deps.Future]]
deps = ["Random"]
uuid = "9fa8497b-333b-5362-9e8d-4d0656e87820"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll", "libdecor_jll", "xkbcommon_jll"]
git-tree-sha1 = "532f9126ad901533af1d4f5c198867227a7bb077"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.4.0+1"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Preferences", "Printf", "Qt6Wayland_jll", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "p7zip_jll"]
git-tree-sha1 = "629693584cef594c3f6f99e76e7a7ad17e60e8d5"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.7"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a8863b69c2a0859f2c2c87ebdc4c6712e88bdf0d"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.7+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "7c82e6a6cd34e9d935e9aa4051b66c6ff3af59ba"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.2+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "d1d712be3164d61d1fb98e7ce9bcbc6cc06b45ed"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.8"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

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

[[deps.InitialValues]]
git-tree-sha1 = "4da0f88e9a39111c2fa3add390ab15f3a44f3ca3"
uuid = "22cec73e-a1b8-11e9-2c92-598750a2cf9c"
version = "0.3.1"

[[deps.IntelOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "14eb2b542e748570b56446f4c50fbfb2306ebc45"
uuid = "1d5cc7b8-4909-519e-a0f8-d0f5ad9712d0"
version = "2024.2.0+0"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "88a101217d7cb38a7b481ccd50d21876e1d1b0e0"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.15.1"
weakdeps = ["Unitful"]

    [deps.Interpolations.extensions]
    InterpolationsUnitfulExt = "Unitful"

[[deps.IntervalSets]]
git-tree-sha1 = "dba9ddf07f77f60450fe5d2e2beb9854d9a49bd0"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.10"
weakdeps = ["Random", "RecipesBase", "Statistics"]

    [deps.IntervalSets.extensions]
    IntervalSetsRandomExt = "Random"
    IntervalSetsRecipesBaseExt = "RecipesBase"
    IntervalSetsStatisticsExt = "Statistics"

[[deps.InverseFunctions]]
git-tree-sha1 = "2787db24f4e03daf859c6509ff87764e4182f7d1"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.16"
weakdeps = ["Dates", "Test"]

    [deps.InverseFunctions.extensions]
    InverseFunctionsDatesExt = "Dates"
    InverseFunctionsTestExt = "Test"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "42d5f897009e7ff2cf88db414a389e5ed1bdd023"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.10.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "c84a835e1a09b289ffcd2271bf2a337bbdda6637"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.3+0"

[[deps.KernelDensity]]
deps = ["Distributions", "DocStringExtensions", "FFTW", "Interpolations", "StatsBase"]
git-tree-sha1 = "7d703202e65efa1369de1279c162b915e245eed1"
uuid = "5ab0869b-81aa-558d-bb23-cbf5423bbe9b"
version = "0.6.9"

[[deps.Kronecker]]
deps = ["LinearAlgebra", "NamedDims", "SparseArrays", "StatsBase"]
git-tree-sha1 = "9253429e28cceae6e823bec9ffde12460d79bb38"
uuid = "2c470bb0-bcc8-11e8-3dad-c9649493f05e"
version = "0.5.5"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "170b660facf5df5de098d866564877e119141cbd"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.2+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "70c5da094887fd2cae843b8db33920bac4b6f07d"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.2+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "5b0d630f3020b82c0775a51d05895852f8506f50"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.4"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LeftChildRightSiblingTrees]]
deps = ["AbstractTrees"]
git-tree-sha1 = "fb6803dafae4a5d62ea5cab204b1e657d9737e7f"
uuid = "1d6d02ad-be62-4b6b-8a6d-2f90e265016e"
version = "0.2.0"

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

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll"]
git-tree-sha1 = "9fd170c4bbfd8b935fdc5f8b7aa33532c991a673"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.11+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "fbb1f2bef882392312feb1ede3615ddc1e9b99ed"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.49.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0c4f9c4f1a50d8f35048fa0532dabbadf702f81e"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.40.1+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "5ee6203157c120d79034c748a2acba45b82b8807"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.40.1+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogDensityProblems]]
deps = ["ArgCheck", "DocStringExtensions", "Random"]
git-tree-sha1 = "f9a11237204bc137617194d79d813069838fcf61"
uuid = "6fdf6af0-433a-55f7-b3ed-c6c6e0b8df7c"
version = "2.1.1"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "a2d09619db4e765091ee5c6ffe8872849de0feea"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.28"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MCMCChains]]
deps = ["AbstractMCMC", "AxisArrays", "Dates", "Distributions", "Formatting", "IteratorInterfaceExtensions", "KernelDensity", "LinearAlgebra", "MCMCDiagnosticTools", "MLJModelInterface", "NaturalSort", "OrderedCollections", "PrettyTables", "Random", "RecipesBase", "Serialization", "Statistics", "StatsBase", "StatsFuns", "TableTraits", "Tables"]
git-tree-sha1 = "c659f7508035a7bdd5102aef2de028ab035f289a"
uuid = "c7f686f2-ff18-58e9-bc7b-31028e88f75d"
version = "5.7.1"

[[deps.MCMCDiagnosticTools]]
deps = ["AbstractFFTs", "DataAPI", "DataStructures", "Distributions", "LinearAlgebra", "MLJModelInterface", "Random", "SpecialFunctions", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "d1737c39191aa26f42a64e320de313f1d1fd74b1"
uuid = "be115224-59cd-429b-ad48-344e309966f0"
version = "0.2.1"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MKL_jll]]
deps = ["Artifacts", "IntelOpenMP_jll", "JLLWrappers", "LazyArtifacts", "Libdl", "oneTBB_jll"]
git-tree-sha1 = "f046ccd0c6db2832a9f639e2c669c6fe867e5f4f"
uuid = "856f044c-d86e-5d09-b602-aeab76dc8ba7"
version = "2024.2.0+0"

[[deps.MLJModelInterface]]
deps = ["Random", "ScientificTypesBase", "StatisticalTraits"]
git-tree-sha1 = "ceaff6618408d0e412619321ae43b33b40c1a733"
uuid = "e80e1ace-859a-464e-9ed9-23947d8ae3ea"
version = "1.11.0"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.MicroCollections]]
deps = ["BangBang", "InitialValues", "Setfield"]
git-tree-sha1 = "629afd7d10dbc6935ec59b32daeb33bc4460a42e"
uuid = "128add7d-3638-4c79-886c-908ea0c25c34"
version = "0.1.4"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.MultivariateStats]]
deps = ["Arpack", "Distributions", "LinearAlgebra", "SparseArrays", "Statistics", "StatsAPI", "StatsBase"]
git-tree-sha1 = "816620e3aac93e5b5359e4fdaf23ca4525b00ddf"
uuid = "6f286f6a-111f-5878-ab1e-185364afe411"
version = "0.10.3"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NamedDims]]
deps = ["LinearAlgebra", "Pkg", "Statistics"]
git-tree-sha1 = "90178dc801073728b8b2d0d8677d10909feb94d8"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "1.2.2"

    [deps.NamedDims.extensions]
    AbstractFFTsExt = "AbstractFFTs"
    ChainRulesCoreExt = "ChainRulesCore"
    CovarianceEstimationExt = "CovarianceEstimation"
    TrackerExt = "Tracker"

    [deps.NamedDims.weakdeps]
    AbstractFFTs = "621f4979-c628-5d54-868e-fcf4e3e8185c"
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    CovarianceEstimation = "587fd27a-f159-11e8-2dae-1979310e6154"
    Requires = "ae029012-a4dd-5104-9daa-d747884805df"
    Tracker = "9f7883ad-71c0-57eb-9f7f-b5c9e6d3789c"

[[deps.NaturalSort]]
git-tree-sha1 = "eda490d06b9f7c00752ee81cfa451efe55521e21"
uuid = "c020b1a1-e9b0-503a-9c33-f039bfc54a85"
version = "1.0.0"

[[deps.NearestNeighbors]]
deps = ["Distances", "StaticArrays"]
git-tree-sha1 = "91a67b4d73842da90b526011fa85c5c4c9343fe0"
uuid = "b8a86587-4115-5ab1-83bc-aa920d37bbce"
version = "0.4.18"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Observables]]
git-tree-sha1 = "7438a59546cf62428fc9d1bc94729146d37a7225"
uuid = "510215fc-4207-5dde-b226-833fc4488ee2"
version = "0.5.5"

[[deps.OffsetArrays]]
git-tree-sha1 = "1a27764e945a152f7ca7efa04de513d473e9542e"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.14.1"
weakdeps = ["Adapt"]

    [deps.OffsetArrays.extensions]
    OffsetArraysAdaptExt = "Adapt"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+4"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a028ee3cb5641cccc4c24e90c36b0a4f7707bdf5"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.14+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "949347156c25054de2db3b166c52ac4728cbad65"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.31"

[[deps.Pango_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "FriBidi_jll", "Glib_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl"]
git-tree-sha1 = "cb5a2ab6763464ae0f19c86c56c63d4a2b0f5bda"
uuid = "36c8627f-9965-5494-a995-c6b170f724f3"
version = "1.52.2+0"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "35621f10a7531bc8fa58f74610b1bfb70a3cfc6b"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.43.4+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "6e55c6841ce3411ccb3457ee52fc48cb698d6fb0"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.2.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PolyaGammaSamplers]]
deps = ["Distributions", "Random", "StatsFuns"]
git-tree-sha1 = "a675f3a525406a8cc4b8d8c024a84e38baddf770"
uuid = "99ff7fc7-8a0f-4729-8284-81f1989d3fc6"
version = "0.1.0"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "66b20dd35966a748321d3b2537c4584cf40387c7"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.3.2"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressLogging]]
deps = ["Logging", "SHA", "UUIDs"]
git-tree-sha1 = "80d919dee55b9c50e8d9e2da5eeafff3fe58b539"
uuid = "33c8b6b6-d38a-422a-b730-caa89a2f386c"
version = "0.1.4"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.PtrArrays]]
git-tree-sha1 = "f011fbb92c4d401059b2212c05c0601b70f8b759"
uuid = "43287f4e-b6f4-7ad1-bb20-aadabca52c3d"
version = "1.2.0"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "492601870742dcd38f233b23c3ec629628c1d724"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.7.1+1"

[[deps.Qt6Declarative_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6ShaderTools_jll"]
git-tree-sha1 = "e5dd466bf2569fe08c91a2cc29c1003f4797ac3b"
uuid = "629bc702-f1f5-5709-abd5-49b8460ea067"
version = "6.7.1+2"

[[deps.Qt6ShaderTools_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll"]
git-tree-sha1 = "1a180aeced866700d4bebc3120ea1451201f16bc"
uuid = "ce943373-25bb-56aa-8eca-768745ed7b5a"
version = "6.7.1+1"

[[deps.Qt6Wayland_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Qt6Base_jll", "Qt6Declarative_jll"]
git-tree-sha1 = "729927532d48cf79f49070341e1d918a65aba6b0"
uuid = "e99dba38-086e-5de3-a5b1-6e4c66e897c3"
version = "6.7.1+1"

[[deps.QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "e237232771fdafbae3db5c31275303e056afaa9f"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.10.1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "1342a47bf3260ee108163042310d26f2be5ec90b"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.5"
weakdeps = ["FixedPointNumbers"]

    [deps.Ratios.extensions]
    RatiosFixedPointNumbersExt = "FixedPointNumbers"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "f65dcb5fa46aee0cf9ed6274ccbd597adc49aa7b"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.7.1"

[[deps.Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e60724fd3beea548353984dc61c943ecddb0e29a"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.4.3+0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.ScientificTypesBase]]
git-tree-sha1 = "a8e18eb383b5ecf1b5e6fc237eb39255044fd92b"
uuid = "30f210dd-8aff-4c5f-94ba-8e64358c1161"
version = "3.0.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.SentinelArrays]]
deps = ["Dates", "Random"]
git-tree-sha1 = "ff11acffdb082493657550959d4feb4b6149e73a"
uuid = "91c51154-3ec4-41a3-a24f-3f23e20d615c"
version = "1.4.5"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Setfield]]
deps = ["ConstructionBase", "Future", "MacroTools", "StaticArraysCore"]
git-tree-sha1 = "e2cc6d8c88613c05e1defb55170bf5ff211fbeac"
uuid = "efcf1570-3423-57d1-acb7-fd33fddbac46"
version = "1.1.1"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.SpecialFunctions]]
deps = ["IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "2f5d4697f21388cbe1ff299430dd169ef97d7e14"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.4.0"
weakdeps = ["ChainRulesCore"]

    [deps.SpecialFunctions.extensions]
    SpecialFunctionsChainRulesCoreExt = "ChainRulesCore"

[[deps.SplittablesBase]]
deps = ["Setfield", "Test"]
git-tree-sha1 = "e08a62abc517eb79667d0a29dc08a3b589516bb5"
uuid = "171d559e-b47b-412a-8079-5efa626c420e"
version = "0.1.15"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "PrecompileTools", "Random", "StaticArraysCore"]
git-tree-sha1 = "eeafab08ae20c62c44c8399ccb9354a04b80db50"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.9.7"
weakdeps = ["ChainRulesCore", "Statistics"]

    [deps.StaticArrays.extensions]
    StaticArraysChainRulesCoreExt = "ChainRulesCore"
    StaticArraysStatisticsExt = "Statistics"

[[deps.StaticArraysCore]]
git-tree-sha1 = "192954ef1208c7019899fbf8049e717f92959682"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.4.3"

[[deps.StatisticalTraits]]
deps = ["ScientificTypesBase"]
git-tree-sha1 = "542d979f6e756f13f862aa00b224f04f9e445f11"
uuid = "64bff920-2084-43da-a3e6-9bb72801c0c9"
version = "3.4.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StatsFuns]]
deps = ["ChainRulesCore", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5950925ff997ed6fb3e985dcce8eb1ba42a0bbe7"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "0.9.18"

[[deps.StatsPlots]]
deps = ["AbstractFFTs", "Clustering", "DataStructures", "Distributions", "Interpolations", "KernelDensity", "LinearAlgebra", "MultivariateStats", "NaNMath", "Observables", "Plots", "RecipesBase", "RecipesPipeline", "Reexport", "StatsBase", "TableOperations", "Tables", "Widgets"]
git-tree-sha1 = "3b1dcbf62e469a67f6733ae493401e53d92ff543"
uuid = "f3b207a7-027a-5e70-b257-86293d7955fd"
version = "0.15.7"

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a04cabe79c5f01f4d723cc6704070ada0b9d46d5"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.3.4"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.TableOperations]]
deps = ["SentinelArrays", "Tables", "Test"]
git-tree-sha1 = "e383c87cf2a1dc41fa30c093b2a19877c83e1bc1"
uuid = "ab02a1b2-a7df-11e8-156e-fb1833f50b87"
version = "1.2.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "OrderedCollections", "TableTraits"]
git-tree-sha1 = "598cd7c1f68d1e205689b1c2fe65a9f85846f297"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.12.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.TerminalLoggers]]
deps = ["LeftChildRightSiblingTrees", "Logging", "Markdown", "Printf", "ProgressLogging", "UUIDs"]
git-tree-sha1 = "f133fab380933d042f6796eda4e130272ba520ca"
uuid = "5d786b92-1e48-4d6f-9151-6b4477ca9bed"
version = "0.1.7"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "96612ac5365777520c3c5396314c8cf7408f436a"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.11.1"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Transducers]]
deps = ["Adapt", "ArgCheck", "BangBang", "Baselet", "CompositionsBase", "ConstructionBase", "DefineSingletons", "Distributed", "InitialValues", "Logging", "Markdown", "MicroCollections", "Requires", "Setfield", "SplittablesBase", "Tables"]
git-tree-sha1 = "3064e780dbb8a9296ebb3af8f440f787bb5332af"
uuid = "28d57a85-8fef-5791-bfe6-a80928e7c999"
version = "0.4.80"

    [deps.Transducers.extensions]
    TransducersBlockArraysExt = "BlockArrays"
    TransducersDataFramesExt = "DataFrames"
    TransducersLazyArraysExt = "LazyArrays"
    TransducersOnlineStatsBaseExt = "OnlineStatsBase"
    TransducersReferenceablesExt = "Referenceables"

    [deps.Transducers.weakdeps]
    BlockArrays = "8e7c35d0-a365-5155-bbbb-fb81a777f24e"
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    LazyArrays = "5078a376-72f3-5289-bfd5-ec5146d43c02"
    OnlineStatsBase = "925886fa-5bf2-5e8e-b522-a9147a512338"
    Referenceables = "42d2dcc6-99eb-4e98-b66c-637b7d73030e"

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

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "d95fe458f26209c66a187b1114df96fd70839efd"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.21.0"
weakdeps = ["ConstructionBase", "InverseFunctions"]

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "975c354fcd5f7e1ddcc1f1a23e6e091d99e99bc8"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.4"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.Widgets]]
deps = ["Colors", "Dates", "Observables", "OrderedCollections"]
git-tree-sha1 = "fcdae142c1cfc7d89de2d11e08721d0f2f86c98a"
uuid = "cc8bc4a8-27d6-5769-a93b-9d913e69aa62"
version = "0.6.6"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "c1a7aa6219628fcd757dede0ca95e245c5cd9511"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "1.0.0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "d9717ce3518dc68a99e6b96300813760d887a01d"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.13.1+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "a54ee957f4c86b526460a720dbc882fa5edcbefc"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.41+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "326b4fea307b0b39892b3e85fa451692eda8d46c"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.1.1+0"

[[deps.Xorg_libSM_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libICE_jll"]
git-tree-sha1 = "3796722887072218eabafb494a13c963209754ce"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.4+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "d2d1a5c49fae4ba39983f63de6afcbea47194e85"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.6+0"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "47e45cd78224c53109495b3e324df0c37bb61fbe"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.11+0"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "bcd466676fef0878338c61e655629fa7bbc69d8e"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.17.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "1827acba325fdcdf1d2647fc8d5301dd9ba43a9d"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.9.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.libdecor_jll]]
deps = ["Artifacts", "Dbus_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pango_jll", "Wayland_jll", "xkbcommon_jll"]
git-tree-sha1 = "9bf7903af251d2050b467f76bdbe57ce541f7f4f"
uuid = "1183f4f0-6f2a-5f1a-908b-139f9cdfea6f"
version = "0.2.2+0"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "490376214c4721cdaca654041f635213c6165cb3"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+2"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.oneTBB_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "7d0ea0f4895ef2f5cb83645fa689e52cb55cf493"
uuid = "1317d2d5-d96f-522e-a858-c73665f53c3e"
version = "2021.12.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─a93d8ec6-7805-11ef-002a-09f4385c2f04
# ╠═8b1cb1c8-921e-471f-a9a4-f02a161b0f11
# ╠═6b9551d9-840d-4afd-b304-d32413fc090d
# ╠═0c3642a4-c7b3-4ded-b490-77c5b73873f1
# ╠═485d8279-ee03-4cef-974e-43c145e94126
# ╠═ea8d8c52-e6aa-44d1-8fdb-7e8f3a4778b5
# ╠═a9602d17-d505-473e-a5f9-90e5ea792437
# ╠═758449bc-5025-4193-8576-bbdeba7c66e3
# ╠═3cc91952-8b1d-4e61-a120-f10c9a5989bd
# ╠═b134cb84-4548-46c4-83e5-caa69c97697d
# ╠═7b12f2aa-a7ab-4dc3-af46-d6260a145a3b
# ╟─4f0acc9a-3c0f-4319-82ce-14f5e6a6a44d
# ╠═cae950e8-75b1-4fd2-8417-7a5d5a9d054d
# ╠═6c23feaa-b1a5-43f8-bd1a-423b18d9eb5f
# ╠═fc521bea-d80c-4381-9df4-3d16aeb74a11
# ╠═c8685063-3ba6-42d0-af83-4b07bc573436
# ╠═84dc174b-1aea-4f55-bb6f-2858a0f6d6b3
# ╠═c7edc678-e100-4bfb-a8e8-55486ba8d225
# ╠═37e82d85-bfeb-482c-a8e4-7ccd1e34431d
# ╠═929f5926-a6be-4db9-b5d7-69154409537f
# ╟─fac11e3a-64f7-4b0c-a3d8-a2c6bfee538b
# ╠═6f1f024f-0047-49b0-86af-c0ef4ceadabf
# ╠═c89252d3-7f4a-471d-8d4f-2db79f9ef4f9
# ╠═2c820be4-7dd5-44fe-8470-32a186f66bae
# ╠═d6cd13fe-0718-43c4-893a-98d7679497d4
# ╠═f2db195d-8ad1-411b-aad5-c2493fe4e431
# ╠═77c2580d-aa1d-436e-a9b2-41aa838ba2ac
# ╠═5bd8ff0d-227b-489d-843c-827708126d7a
# ╠═14786bb2-28aa-498e-b7f1-eadcf4718a0e
# ╟─cdc7fbd2-8315-4c50-832f-0e1e353f1172
# ╟─5a502f85-2624-4ba9-abea-7428f1989816
# ╠═1000ab97-3d0b-4f3f-984a-c7342d0e79d8
# ╠═fbed8560-7adc-4ad4-a14f-4dbda25ad67e
# ╠═970a24bf-8ce5-4e9e-89a9-fc7a6c24d36b
# ╠═0f6abf0d-aab8-4b40-806f-978f9cbb2f37
# ╠═9c437f8f-d1bc-4218-9f83-5a438c8c22b6
# ╟─a975df83-7e0f-45eb-bb1b-c2c673d5d6e9
# ╠═de279a3f-0c73-4542-8aaa-9d397d2c8196
# ╠═8ba1f85f-7498-4b3e-82e2-36bf29998265
# ╠═718c896b-5746-4ad9-9e0b-fa7af2dc56cd
# ╠═c7c96185-af24-4c38-b7af-c1897acd99bd
# ╟─888f5005-a524-42e6-b8a2-879418e4248a
# ╠═89f381ed-f0fb-462e-b014-805c1575e139
# ╠═1f9287ca-cb85-4e7a-8e5f-32a343f94e2c
# ╠═9a6816dd-fcd4-4fe8-a8e0-f17038976077
# ╠═94583a08-4dfb-4579-b8a2-539833c6645d
# ╠═26771447-5b21-4733-95d5-38609d0239d8
# ╠═07522aa2-ce40-4a90-877c-ff1b5fadd789
# ╠═c5c72069-ab2b-4f79-a80b-d256441a8a6c
# ╟─6c19bce2-70f2-43d9-9fdf-8cf1bb505ac1
# ╠═91235ba9-0fb9-4452-86a3-d4a88db41b7e
# ╠═ad33ce91-0f89-4520-89fb-3f06645d4df3
# ╠═2c4f365d-321b-4588-a49d-5a29dc4432b0
# ╠═ea14de6e-15b0-4c3f-af7b-cdf57de0bf2e
# ╟─30f39851-d610-4faa-adf8-1b81a42c09a5
# ╠═e1c7895c-a770-4e97-bb21-04a34baaa14c
# ╠═4626846f-4c5a-4171-bcec-9d64bcc07bac
# ╠═aaa4e2de-90c9-43c3-b36d-8eb23842419f
# ╠═8a383a20-8a9b-428b-8e11-d8e3229bd3ff
# ╠═762a1217-6892-4b1e-9ee4-ed6de7f63f7e
# ╠═335b4fc4-a079-4107-8726-2fa802468041
# ╟─5eeb8329-0c72-49b7-9179-07df457861fc
# ╠═e5c47428-e543-4836-8f93-57ea61aec05d
# ╠═39ea5fbf-3d06-4050-b8b1-12c9eddff127
# ╠═0c3be9d9-89d1-4184-b80a-adbe8c86d15c
# ╠═563c6957-d66c-4b74-8442-ad0c74aeff9f
# ╟─3fe057b5-7dda-4843-8f58-0cde265cfc80
# ╟─d4b00edb-7a00-480f-a8dd-280daf5087e1
# ╠═6a8d0b9b-c11d-4f79-b2aa-80d47b3b8094
# ╠═1af30d48-3553-4fe1-8a2a-b5c598740b71
# ╟─039c4bba-e87f-430c-bbe1-8e6c2d98df02
# ╠═38b93062-57e4-486b-af9c-aa4a73080f01
# ╠═49b24a83-2385-4f9d-a607-fd685a953d89
# ╠═b9eb3646-1f58-4a14-8615-c5850a31f65d
# ╠═2ba37c20-f1e1-4875-b9ec-6ba0647f5b1b
# ╠═a59e9814-60e9-4c91-ac22-3883eb154c9d
# ╠═8259b697-54e3-448a-9332-4ad9e7c7fc45
# ╠═5d1d2385-bba0-4ac4-8336-938078c1a7e8
# ╠═b98d9dce-6c4f-4a6f-b4d4-19b42bf96c08
# ╠═14a33fb9-0733-4db8-90b9-89c370458b98
# ╠═5cdb09bd-8e35-4b46-a14f-b32de5ac0429
# ╠═76767de6-85d2-4af7-9c54-3c4f0070002d
# ╠═e4564810-15fc-4820-9da1-b0818bd6be0e
# ╠═4c88ac80-d6a4-4192-a4e8-e4bc23fecc59
# ╠═708eb634-9e44-4164-9afd-33225c8a9013
# ╠═fcba6869-3224-46e3-9b3b-b754e7671499
# ╠═44c6092d-9a02-4718-8e44-db407b5f124e
# ╠═8c3cee3e-e283-49d2-9425-321c128d528a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

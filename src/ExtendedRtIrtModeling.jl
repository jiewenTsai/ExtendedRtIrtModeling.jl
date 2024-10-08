
"""
# An Extended RtIrt Modeling

In this note, I primarily collect the structs and functions used in QrSA model.

"""

#using ProgressMeter

"""

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


# ╔═╡ 485d8279-ee03-4cef-974e-43c145e94126
"""
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
"""
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
"""
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
"""
## 1.2 _draw_-Functions
"""

# ╔═╡ 5a502f85-2624-4ba9-abea-7428f1989816
" ### 1.2.1 2PL IRT"

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
"""
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
"""
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
"""

## ⭐ 打包區

"""

# ╔═╡ 30f39851-d610-4faa-adf8-1b81a42c09a5
"""
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

### A Pluto.jl notebook ###
# v0.20.3

#using Markdown
#using InteractiveUtils

# ╔═╡ a40a31e9-cc13-4d36-9186-212e8a562ddf
using PlutoUI,
	ProgressMeter,
	PrettyTables

# ╔═╡ 076769b0-9f55-11ef-28c9-4fbd31a7b602
md"""
# GibbsRtIrt


- GibbsRtIrtNull
- GibbsRtIrt
"""

# ╔═╡ d99a9ef6-8ab0-42a2-9e5c-8c6f4c591e9e
TableOfContents()

# ╔═╡ a1da4c97-859e-4d3e-a131-a5f50dd3e102
md"""
## Structs
"""

# ╔═╡ 709be8fe-d484-4371-97d8-f41e02113b91
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

# ╔═╡ d7f4e8c6-2593-4b80-b7e4-67300135e2d4
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

# ╔═╡ f907e709-5933-4a19-a3fe-7f3abee405c3
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

# ╔═╡ 1653e19e-b0f4-4ebc-aa7b-44678648e36f
abstract type GibbsRtIrt2 end

# ╔═╡ 6e9a1916-137c-4d6d-90ef-5b6e3d0d8a24
"""
"""
mutable struct GibbsRtIrt <: GibbsRtIrt2
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
			ζ=randn(self.Cond.nSubj),
			λ=zeros(self.Cond.nItem),
			σ²t=ones(self.Cond.nItem),
			β=randn(self.Cond.nFeat+1,2), 
			Σp=I(2) )
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

# ╔═╡ a213c2de-370d-4580-8dd2-9b1c959d6393
"""
"""
mutable struct GibbsRtIrtNull <: GibbsRtIrt2
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
			ζ=randn(self.Cond.nSubj),
			λ=zeros(self.Cond.nItem),
			σ²t=ones(self.Cond.nItem),
			#β=randn(self.Cond.nFeat,2), 
			Σp=I(2) )
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


# ╔═╡ 53ff336c-0460-4473-a426-e856783a7fd9
md"""
## Functions
"""

# ╔═╡ bbce8f4c-82e0-437a-850c-c2e8ff63f954
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

# ╔═╡ 03d68d77-ff66-4e85-8e46-fa9a1d42cdff
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

# ╔═╡ 1c157bce-8ef6-47f2-a8b5-37997cffb87a
"""
"""
function getLogLikelihoodRtIrt(Cond,Data; P=Para)
	x = [ones(Cond.nSubj) Data.X]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.λ' .- P.ζ
	η = [P.θ P.ζ]
    μη = x * reshape(P.β, Cond.nFeat+1, 2)
	Ση = reshape(P.Σp, 2,2) 
	logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, sqrt.(P.σ²t')), Data.logT)) sum([logpdf(MvNormal(μη[i,:], Ση), η[i,:]) for i in 1:Cond.nSubj])]

    return sum(logPdf)
end

# ╔═╡ 4d78cd12-de4b-44c9-8477-cdcd88a368ac
"""
	sample(MCMC::GibbsRtIrt)
"""
function sample!(MCMC::GibbsRtIrt; intercept=false, itemtype::Union{String} = "2pl", cov2one=true)

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	

		## structural
		Para.β = drawSubjCoefficients(Cond, Data, Para )
		if intercept == false
			Para.β[1,:] = [0. 0.]
		end
		Para.Σp = drawSubjCovariance(Cond, Data, Para, cov2one)
		

		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.b = drawItemDifficulty(Data,Para)
		Para.a = drawItemDiscrimination(Data, Para)
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		Para.θ = drawSubjAbility(Cond, Data, Para)

		
		
		## response time part
	    Para.λ = drawItemIntensity(Cond,Data,Para)
	    Para.σ²t = drawItemTimeResidual(Cond,Data,Para)
		Para.ζ = drawSubjSpeed(Cond,Data,Para)




		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.ζ; Para.λ; Para.σ²t]
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
		ζ = mean(Post.rt[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		λ = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		σ²t = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec

	)

	return MCMC
end

# ╔═╡ d9eaeeb6-aa70-4142-b1ee-2520d4391b45
"""
"""
function getLogLikelihoodRtIrtNull(Cond,Data; P=Para)
	#x = [ones(Cond.nSubj) Data.X]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.λ' .- P.ζ
	η = [P.θ P.ζ]
    μη = zeros(Cond.nSubj, 2) #x * reshape(P.β, Cond.nFeat+1, 2)
	Ση = reshape(P.Σp, 2,2)
	logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, sqrt.(P.σ²t')), Data.logT)) sum([logpdf(MvNormal(μη[i,:], Ση), η[i,:]) for i in 1:Cond.nSubj])]

    return sum(logPdf)
end

# ╔═╡ 97a5de5d-8593-4d9b-bdcf-26cf05d97047
"""
	sample(MCMC::GibbsRtIrtNull)
"""
function sample!(MCMC::GibbsRtIrtNull; itemtype::Union{String} = "2pl", cov2one=true)

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	
		## structural
		Para.β = zeros(Cond.nFeat+1, 2)
		Para.Σp = drawSubjCovarianceNull(Cond, Data, Para, cov2one)
		
		## irt part
		Para.ω = drawRaPgRandomVariable(Para)
		Para.b = drawItemDifficulty(Data,Para)
		Para.a = drawItemDiscrimination(Data, Para)	
		if itemtype == "1pl"
			Para.a = ones(Cond.nItem)
		end
		
		Para.θ = drawSubjAbilityNull(Cond, Data, Para)
		
		## response time part
	    Para.λ = drawItemIntensity(Cond,Data,Para)
	    Para.σ²t = drawItemTimeResidual(Cond,Data,Para)
		Para.ζ = drawSubjSpeedNull(Cond,Data,Para)

		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.ζ; Para.λ; Para.σ²t]
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
		ζ = mean(Post.rt[(Cond.nBurnin+1):end, 1:(Cond.nSubj), :], dims=(1,3)) |> vec,
		λ = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+1):(Cond.nSubj+Cond.nItem), :], dims=(1,3)) |> vec,
		σ²t = mean(Post.rt[(Cond.nBurnin+1):end, (Cond.nSubj+Cond.nItem+1):end, :], dims=(1,3)) |> vec

	)

	return MCMC
end

# ╔═╡ fea66a6d-6f90-45e7-b15d-f3c0858257a1
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

# ╔═╡ d43abe71-20ee-4923-aca3-2388c837bbef
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

# ╔═╡ e422ca61-6931-4818-bd4a-064fcb232eda
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


# ╔═╡ fe805a42-b5da-47bb-a6fe-170426b1099a
"""
	coef(MCMC)
"""
function coef(MCMC::GibbsRtIrt2)
    h1 = Highlighter(
        f = (data, i, j) -> (data[i, j] isa AbstractFloat && data[i, j] <= 0.),
        crayon = crayon"red"
    )
    dfItem =  DataFrame(
        Item = collect(1:MCMC.Cond.nItem), 
        a = MCMC.Post.mean.a, 
        b = MCMC.Post.mean.b,
        λ = MCMC.Post.mean.λ,
        σ²t = MCMC.Post.mean.σ²t
    )
    CovIndex = ["θ", "ζ"]
    dfCov = DataFrame(
        [CovIndex reshape(MCMC.Post.mean.Σp, 2,2)],
        :auto
    )

    βIndex = ["β$i" for i in 0:MCMC.Cond.nFeat]
    dfBeta = DataFrame(
        [βIndex reshape(MCMC.Post.mean.β, MCMC.Cond.nFeat+1,2)],
        :auto
    )
    DIC = getDic(MCMC)

    ## display
	println("=="^30)
	println(">> Model: $(typeof(MCMC)).")
	println(">> This analysis involved $(MCMC.Cond.nSubj) subjects, $(MCMC.Cond.nItem) items, and $(MCMC.Cond.nFeat) features.")
	println(">> $(MCMC.Cond.nChain) chains of $(MCMC.Cond.nIter) iterations each were run,")
	if MCMC.Cond.nThin > 1
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, and every $(MCMC.Cond.nThin) iterations kept,")
		println(">> totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	else
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	end
	
	if typeof(MCMC) == "GibbsRtIrtQr"
		println(">> RT Quantile: $(MCMC.Cond.qRt)")
	end
	println("=="^30)

	
    println("1) Item Parameters.")
    pretty_table(dfItem, highlighters=h1,formatters = ft_printf("%5.3f", 2:5))

    println("2) Covariance of Person Parameters.")
    pretty_table(dfCov, header=["Coef", "θ", "ζ"], highlighters=h1,formatters = ft_printf("%5.3f"))

    println("3) Regression Coefficients.")
    pretty_table(dfBeta, header=["Coef", "θ", "ζ"], highlighters=h1, formatters = ft_printf("%5.3f") )

    println("4) Criterion.")
    pretty_table(
        [DIC.DIC-DIC.pD DIC.DIC], header=["Deviance", "DIC"], 
        highlighters=h1,
        formatters = ft_printf("%5.3f")
    )

end

# ╔═╡ a532a975-6975-4125-a689-3d31502bbfba
"""
	coef(MCMC::GibbsMlIrt)
"""
function coef(MCMC::GibbsMlIrt)
    h1 = Highlighter(
        f = (data, i, j) -> (data[i, j] isa AbstractFloat && data[i, j] <= 0.),
        crayon = crayon"red"
    )
    dfItem =  DataFrame(
        Item = collect(1:MCMC.Cond.nItem), 
        a = MCMC.Post.mean.a, 
        b = MCMC.Post.mean.b,
    )
    
    CovIndex = ["θ"]
    dfCov = DataFrame([CovIndex 1.], :auto)
    

    βIndex = ["β$i" for i in 0:MCMC.Cond.nFeat]
    dfBeta = DataFrame(
        [βIndex reshape(MCMC.Post.mean.β, MCMC.Cond.nFeat+1,1)],
        :auto
    )

    DIC = getDic(MCMC)

    ## display
	println("=="^30)
	println(">> Model: $(typeof(MCMC)).")
	println(">> This analysis involved $(MCMC.Cond.nSubj) subjects, $(MCMC.Cond.nItem) items, and $(MCMC.Cond.nFeat) features.")
	println(">> $(MCMC.Cond.nChain) chains of $(MCMC.Cond.nIter) iterations each were run,")
	if MCMC.Cond.nThin > 1
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, and every $(MCMC.Cond.nThin) iterations kept,")
		println(">> totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	else
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	end
	
	if typeof(MCMC) == "GibbsRtIrtQr"
		println(">> RT Quantile: $(MCMC.Cond.qRt)")
	end
	println("=="^30)

	
    println("1) Item Parameters.")
    pretty_table(dfItem, highlighters=h1, formatters = ft_printf("%5.3f", 2:3))

    println("2) Covariance of Person Parameters.")
    pretty_table(dfCov, header=["Coef", "θ"], highlighters=h1, formatters = ft_printf("%5.3f"))

    println("3) Regression Coefficients.")
    pretty_table(dfBeta, header=["Coef", "θ"], highlighters=h1, formatters = ft_printf("%5.3f"))

    println("4) Criterion.")
    pretty_table(
        [DIC.DIC-DIC.pD DIC.DIC], header=["Deviance", "DIC"], 
        highlighters=h1,
        formatters = ft_printf("%5.3f")
    )

end

# ╔═╡ ec8da7d9-8f20-4f92-986e-1f8690be2039
"""
	precis(MCMC)
"""
function precis(MCMC::GibbsRtIrt2)

    mcmcRa = MCMC.Post.ra[(MCMC.Cond.nBurnin+1):end, (MCMC.Cond.nSubj+1):end, :]
    chainsMcmcRa = Chains(mcmcRa, [["a$i" for i in 1:MCMC.Cond.nItem]; ["b$i" for i in 1:MCMC.Cond.nItem]])

    mcmcRt = MCMC.Post.rt[(MCMC.Cond.nBurnin+1):end, (MCMC.Cond.nSubj+1):end, :]
    chainsMcmcRt = Chains(mcmcRt, [["λ$i" for i in 1:MCMC.Cond.nItem]; ["σ²t$i" for i in 1:MCMC.Cond.nItem]])

    mcmcQr = MCMC.Post.qr[(MCMC.Cond.nBurnin+1):end, :, :]
    chainsMcmcQr = Chains(mcmcQr, [vec(["β[$i,$j]" for i in 0:(MCMC.Cond.nFeat), j in 1:2]); ["Σ[1,1]","Σ[1,2]","Σ[2,1]","Σ[2,2]"]])

    ## display
	println("=="^30)
	println(">> Model: $(typeof(MCMC)).")
	println(">> This analysis involved $(MCMC.Cond.nSubj) subjects, $(MCMC.Cond.nItem) items, and $(MCMC.Cond.nFeat) features.")
	println(">> $(MCMC.Cond.nChain) chains of $(MCMC.Cond.nIter) iterations each were run,")
	if MCMC.Cond.nThin > 1
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, and every $(MCMC.Cond.nThin) iterations kept,")
		println(">> totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	else
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	end
	
	if typeof(MCMC) == "GibbsRtIrtQr"
		println(">> RT Quantile: $(MCMC.Cond.qRt)")
	end
	println("=="^30)

	
    println("1) Item Response Model.")
    getPrecisTable(chainsMcmcRa)
    println("2) Response Time Model.")
    getPrecisTable(chainsMcmcRt)
    println("3) Structural Model.")
    getPrecisTable(chainsMcmcQr)

end

# ╔═╡ 51ce48f3-bdef-4642-a814-edd0f2c074cd
"""
	precis(MCMC::GibbsMlIrt)
"""
function precis(MCMC::GibbsMlIrt)

    mcmcRa = MCMC.Post.ra[(MCMC.Cond.nBurnin+1):end, (MCMC.Cond.nSubj+1):end, :]
    chainsMcmcRa = Chains(mcmcRa, [["a$i" for i in 1:MCMC.Cond.nItem]; ["b$i" for i in 1:MCMC.Cond.nItem]])

    ## display
	println("=="^30)
	println(">> Model: $(typeof(MCMC)).")
	println(">> This analysis involved $(MCMC.Cond.nSubj) subjects, $(MCMC.Cond.nItem) items, and $(MCMC.Cond.nFeat) features.")
	println(">> $(MCMC.Cond.nChain) chains of $(MCMC.Cond.nIter) iterations each were run,")
	if MCMC.Cond.nThin > 1
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, and every $(MCMC.Cond.nThin) iterations kept,")
		println(">> totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	else
		println(">> with the first $(MCMC.Cond.nBurnin) discarded, totaling $(Int(MCMC.Cond.nChain * (MCMC.Cond.nIter - MCMC.Cond.nBurnin) / MCMC.Cond.nThin)) saved iterations.")
	end
	
	if typeof(MCMC) == "GibbsRtIrtQr"
		println(">> RT Quantile: $(MCMC.Cond.qRt)")
	end
	println("=="^30)

	
    println("1) Item Response Model.")
    getPrecisTable(chainsMcmcRa)    

end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PrettyTables = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"

[compat]
PlutoUI = "~0.7.60"
PrettyTables = "~2.4.0"
ProgressMeter = "~1.10.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.3"
manifest_format = "2.0"
project_hash = "6a91486a11a55c4e1a004aeab89a1fa6ef02e9a4"

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

[[deps.Crayons]]
git-tree-sha1 = "249fe38abf76d48563e2f4556bebd215aa317e15"
uuid = "a8cc5b0e-0ffa-5ad4-8c14-923d3ee1735f"
version = "4.1.1"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

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

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

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

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

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

[[deps.PrettyTables]]
deps = ["Crayons", "LaTeXStrings", "Markdown", "PrecompileTools", "Printf", "Reexport", "StringManipulation", "Tables"]
git-tree-sha1 = "1101cd475833706e4d0e7b122218257178f48f34"
uuid = "08abe8d2-0d0c-5749-adfa-8a2ac140af0d"
version = "2.4.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

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

[[deps.StringManipulation]]
deps = ["PrecompileTools"]
git-tree-sha1 = "a6b1675a536c5ad1a60e5a5153e1fee12eb146e3"
uuid = "892a3eda-7b42-436c-8928-eab12a02cf0e"
version = "0.4.0"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

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
# ╟─076769b0-9f55-11ef-28c9-4fbd31a7b602
# ╠═a40a31e9-cc13-4d36-9186-212e8a562ddf
# ╠═d99a9ef6-8ab0-42a2-9e5c-8c6f4c591e9e
# ╠═a1da4c97-859e-4d3e-a131-a5f50dd3e102
# ╠═709be8fe-d484-4371-97d8-f41e02113b91
# ╠═d7f4e8c6-2593-4b80-b7e4-67300135e2d4
# ╠═f907e709-5933-4a19-a3fe-7f3abee405c3
# ╠═1653e19e-b0f4-4ebc-aa7b-44678648e36f
# ╠═6e9a1916-137c-4d6d-90ef-5b6e3d0d8a24
# ╠═a213c2de-370d-4580-8dd2-9b1c959d6393
# ╟─03d68d77-ff66-4e85-8e46-fa9a1d42cdff
# ╠═4d78cd12-de4b-44c9-8477-cdcd88a368ac
# ╠═97a5de5d-8593-4d9b-bdcf-26cf05d97047
# ╟─53ff336c-0460-4473-a426-e856783a7fd9
# ╠═bbce8f4c-82e0-437a-850c-c2e8ff63f954
# ╠═1c157bce-8ef6-47f2-a8b5-37997cffb87a
# ╠═d9eaeeb6-aa70-4142-b1ee-2520d4391b45
# ╠═fea66a6d-6f90-45e7-b15d-f3c0858257a1
# ╠═d43abe71-20ee-4923-aca3-2388c837bbef
# ╠═e422ca61-6931-4818-bd4a-064fcb232eda
# ╠═d97fb13b-35e2-40fd-aedf-dd72db2d0e90
# ╠═fe805a42-b5da-47bb-a6fe-170426b1099a
# ╠═a532a975-6975-4125-a689-3d31502bbfba
# ╠═ec8da7d9-8f20-4f92-986e-1f8690be2039
# ╠═51ce48f3-bdef-4642-a814-edd0f2c074cd
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

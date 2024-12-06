### A Pluto.jl notebook ###
# v0.20.3

#using Markdown
#using InteractiveUtils

# ╔═╡ 1d4a9f51-3db8-4a78-936d-a7f2879cc0a0
using PlutoUI

# ╔═╡ 3dd409f0-ab34-11ef-22d1-3b72e5416118
md"""
# Simulation Study Tools


## prepare real data.
"""

# ╔═╡ 23484455-73cd-4c7f-a13d-f58fb0bc3051
TableOfContents()

# ╔═╡ bc0d18af-ab68-443d-8e98-ec5fbedc2e12
mutable struct OutputMetrics
	Key
    Rmse 
    Bias 
    Corr 
    Dic 
	Diag
    function OutputMetrics(;
		Key = [],
        Rmse = [],
        Bias = [],
        Corr = [],
        Dic = [],
		Diag = []
    )
        return new(Key, Rmse, Bias, Corr, Dic, Diag)
    end
end

# ╔═╡ b57756d7-b28a-4fba-be7d-d1b9f196fb20
begin
	getRmse(a,b) = mean(sqrt(mean((a - b).^2)))
	getBias(a,b) = mean(mean(a - b))
end

# ╔═╡ d66e453f-de50-42a5-a821-508e23b44b10
"""
begin
	function getPropertyPost(MCMC, name)
		Base.getproperty(MCMC, :Post) |> 
		    x -> Base.getproperty(x, :mean) |>
		    x -> Base.getproperty(x, name) 
	end
	function getPropertyTrue(MCMC, name)
		Base.getproperty(MCMC, :truePara) |> 
		    x -> Base.getproperty(x, name) 
	end
end
"""

# ╔═╡ 3468403a-3ae8-4635-9fd9-66538be96315
"""
begin
	getRmse(MCMC, name) = mean(sqrt(mean((getPropertyPost(MCMC, name) .- getPropertyTrue(MCMC, name)).^2)))
	getBias(MCMC, name) = mean(mean(getPropertyPost(MCMC, name) .- getPropertyTrue(MCMC, name)))
	getCorr(MCMC, name) = cor(getPropertyPost(MCMC, name), getPropertyTrue(MCMC, name))
end
"""

# ╔═╡ 74c8ac72-0701-4a5a-8409-0d93e6022503
"""
"""
function setTrueParaRtIrt(Cond;
	a = [],
	b = [],
	β = [],
	λ = [],
	σ²t = [],
	trueStdRa=1.,
	trueStdRt=1.,
	trueCorr=0.
	
)
	truePara = InputPara()
	truePara.a = rand(Truncated.(Normal(1., (0.2)),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., (0.5)), Cond.nItem)
	truePara.λ = rand(Truncated(Normal(4., (0.2) ), 0, Inf), Cond.nItem)
	truePara.σ²t = rand(LogNormal(log(0.3), (0.2)), Cond.nItem)
	trueStd = Diagonal([trueStdRa, trueStdRt])
    trueCor = [1. trueCorr; trueCorr 1.]
	truePara.Σp = trueStd * trueCor * trueStd
	truePara.β = rand(MvNormal(zeros(2), I(2)), Cond.nFeat)'
	return truePara 
end

# ╔═╡ 201f594a-cc89-4236-957f-ed88953b9ebd
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
	truePara.a = rand(Truncated.(Normal(1., (0.2)),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., (0.5)), Cond.nItem)
	truePara.β = rand(MvNormal(zeros(1), I(1)), Cond.nFeat)'
	return truePara 
end

# ╔═╡ 595dbac9-1bb6-4975-928c-5a53f28eae56
"""
"""
function setDataRtIrtNull(Cond, truePara)

	## structure

	## ra and rt
	noise = rand(MvNormal(zeros(2), truePara.Σp), Cond.nSubj)
    trueSubj = noise'
	truePara.θ = trueSubj[:,1]
	truePara.ζ = trueSubj[:,2]
	#truePara.Σp = cov(trueSubj')


	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    μt =  truePara.λ' .- truePara.ζ 
    trueT = rand.(Truncated.(LogNormal.(μt, sqrt.(truePara.σ²t')), 0, Inf))

	## collect
	trueData = InputData(Y=trueY, T=trueT)


	return trueData

end

# ╔═╡ 838de637-d131-480b-99d8-ed55433f41de
"""
"""
function setDataRtIrt(Cond, truePara)

	## structure
	trueX = Array{Float64}(undef, Cond.nSubj, Cond.nFeat)
	trueX[:,1:end] = rand(Normal(0, 1.), Cond.nSubj, (Cond.nFeat))

	## ra and rt
    trueMean = trueX * truePara.β
	#noise = cholesky(Symmetric(truePara.Σp)).L * randn(2,Cond.nSubj)
	noise = rand(MvNormal(zeros(2), truePara.Σp), Cond.nSubj)
    trueSubj = trueMean .+ noise'
	truePara.θ = trueSubj[:,1]
	truePara.ζ = trueSubj[:,2]
	
	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    μt =  truePara.λ' .- truePara.ζ
    trueT = rand.(Truncated.(LogNormal.(μt, sqrt.(truePara.σ²t')), 0, Inf))

	## collect
	trueData = InputData(Y=trueY, X=trueX, T=trueT)


	return trueData

end

# ╔═╡ e3e68b65-9e25-49c7-8a8f-07b34dfa1b0f
md"""
### Cross and Latent
"""

# ╔═╡ 732cc1f5-549e-492f-a37d-86ed70cdb7b8
"""
"""
function setTrueParaRtIrtCross(Cond;
	a = [],
	b = [],
	ρ = [],
	λ = [],
	σ²t = [],
	trueStdRa=1.,
	trueStdRt=1.
	
)
	truePara = InputPara()
	truePara.a = rand(Truncated.(Normal(1., (0.2)),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., (0.5)), Cond.nItem)
	truePara.λ = rand(Truncated(Normal(3., (0.5) ), 0, Inf), Cond.nItem)
	truePara.σ²t = rand(LogNormal(log(0.3), (0.2)), Cond.nItem)
	trueStd = Diagonal([trueStdRa, trueStdRt])
    trueCor = I(2)
	truePara.Σp = trueStd * trueCor * trueStd
	truePara.ρ = rand(Normal( 0.2, (0.1)), Cond.nItem)
	#truePara.ρ = rand(Normal( 0.3, (0.1)), Cond.nItem)
	return truePara 
end

# ╔═╡ f8e59e9a-16d6-4e4e-b630-f81f593d9f01

"""
	setDataRtIrtCross(Cond, truePara; type)

### Arguments.

	- type: "norm", "tail", "skew"
"""
function setDataRtIrtCross(Cond, truePara; type="norm")

	## structure

	## ra and rt
	noise = rand(MvNormal(zeros(2), truePara.Σp), Cond.nSubj)
    trueSubj = noise'
	truePara.θ = trueSubj[:,1]
	truePara.ζ = trueSubj[:,2]

	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    μt =  truePara.λ' .- truePara.ζ .- truePara.θ  * (truePara.ρ')

    if type == "norm"
		#trueT = μt .+ rand(Normal(0, 0.5), Cond.nSubj, Cond.nItem)
        trueT = rand.(Normal.(μt, sqrt.(truePara.σ²t')))
    elseif type == "tail"
		#trueT = μt .+ rand.(Normal.(0, sqrt.(truePara.σ²t'))) 
		#trueT = μt .+ sqrt.(truePara.σ²t') .* rand.(TDist(3))	
		trueT = μt .+ rand(TDist(3), Cond.nSubj, Cond.nItem)
		#trueT =  rand.(Truncated.(Cauchy.(μt, 1),-6,6))
		# 對於尾部分配
		#σ_adjusted = sqrt.(truePara.σ²t' .* π/2)  # 調整Cauchy尺度
		#trueT = μt .+ rand.(Truncated.(Cauchy.(0, σ_adjusted),-6,6))
        
		#trueT = μt .+ rand.(Truncated.(Cauchy.(0, sqrt.(truePara.σ²t')),-6,6))
        #trueT = μt .+ abs.(μt).*sqrt.(truePara.σ²t').*randn()
        #trueT = μt .+ 0.5*abs.(μt).*randn()

		# 應該改為：
		#ν = 4
		#scaling = sqrt.((ν-2)/ν)  # 調整變異數
		#trueT = μt .+ sqrt.(truePara.σ²t') .* scaling .* rand.(TDist(ν))


    elseif type == "skew"
		trueT = μt .+ rand.(Gamma.(1, 0.5)) .- 1
		#trueT = rand.(LogNormal.( log.(μt), sqrt.(truePara.σ²t')))
        #trueT = μt .+ rand.(Gamma.(sqrt.(truePara.σ²t'), 1)) .- 0.5
		#trueT = μt .+ rand(Gamma(0.5, 1), Cond.nSubj, Cond.nItem) .- 0.5
		#trueT = μt .+ rand(Beta(0.5,0.5), Cond.nItem)'

		# 應該改為：
		#σ_log = 0.8  # 控制偏度大小
		#noise = rand.(LogNormal.(0, sqrt.(truePara.σ²t')), Cond.nSubj)
		#trueT = μt .+ noise
		
    end
    #trueT = clamp.(trueT, 0,6)
    trueT = exp.(trueT)
	#trueT = clamp.(trueT, 0, exp(6) )

	## collect
	trueData = InputData(Y=trueY, T=trueT)


	return trueData

end

# ╔═╡ 54a39e13-807d-49ba-9cd0-710c7e500eff
"""
"""
function setTrueParaRtIrtLatent(Cond;
	a = [],
	b = [],
	β = [],
	λ = [],
	σ²t = [],
	trueStdRa=1.,
	trueStdRt=1.
	
)
	truePara = InputPara()
	truePara.a = rand(Truncated.(Normal(1., (0.2)),0,Inf), Cond.nItem)
	truePara.b = rand(Normal(0., (0.5)), Cond.nItem)
	truePara.λ = rand(Truncated(Normal(3., (0.2) ), 0, Inf), Cond.nItem)
	truePara.σ²t = rand(Truncated.(Normal(0.3, (0.5)), 0, Inf), Cond.nItem)
	trueStd = Diagonal([trueStdRa, trueStdRt])
    trueCor = I(2)
	truePara.Σp = trueStd * trueCor * trueStd
    trueρ = rand(Truncated(Normal(0,0.5),-1,1))
	truePara.β = [rand( Normal( 0, (0.5)), Cond.nFeat); trueρ]
	return truePara 
end


# ╔═╡ 7fa851f0-f16d-43d1-af39-a992f5764c88

"""
	setDataRtIrtLatent(Cond, truePara; type)

### Arguments.

	- type: "norm", "tail", "skew"
"""
function setDataRtIrtLatent(Cond, truePara; type="norm")

	## structure

	## ra and rt
	truePara.θ = randn(Cond.nSubj)

    trueX = rand(Normal(0, 1.), Cond.nSubj, Cond.nFeat)
    x = [trueX truePara.θ]


    ## rt
    μt = x * truePara.β #zeros(Cond.nSubj) # 

    if type == "norm"
        #truePara.ζ = μt .+ randn(Cond.nSubj)
		truePara.ζ =  rand.(Normal.(μt, 0.5))
    elseif type == "tail"

		truePara.ζ = μt .+ rand(TDist(3), Cond.nSubj)
        #truePara.ζ = μt .+ rand(Truncated(Cauchy(0, 1),-6,6),Cond.nSubj)
        #trueT = μt .+ abs.(μt).*sqrt.(truePara.σ²t').*randn()
        #trueT = μt .+ 0.5*abs.(μt).*randn()
    elseif type == "skew"
        truePara.ζ = μt .+ rand(Gamma(1, 0.5), Cond.nSubj) .- 1
        #truePara.ζ = rand.(LogNormal.(μt, 0.25)) #.- exp.(μt .+ 0.25^2/2)
    end


	## irt
    truePr =  truePara.a' .* (truePara.θ .- truePara.b') 
    trueY = rand.(BernoulliLogit.(truePr))

	## rt
    trueμt =  truePara.λ' .- truePara.ζ
    trueT = rand.(LogNormal.(trueμt, sqrt.(truePara.σ²t')))

	## collect
	trueData = InputData(Y=trueY, T=trueT, X=trueX)

	return trueData

end


# ╔═╡ 7b1041d1-bb69-412a-80e8-c0ab731545f4
"""
"""
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

# ╔═╡ cf125fa8-18c0-4f44-80fd-4447c45a5876
"""
"""
function testingDict(nSubj::Int, nItem::Int , nFeat::Int)
	data = Dict(
		"Y" => rand(nSubj, nItem),
		"" => rand(nSubj, nFeat),
		"T" => randn(nSubj, nItem)
	
	)
	return data
end


# ╔═╡ d792cfdb-503f-4986-a62d-02338db36ac4

"""
	comparePara(Mcmc; name=:a)
"""
function comparePara(Mcmc; name=:a)

    if name==:β
        True = getfield(Mcmc.truePara, Symbol(name))
        Esti = vec(getfield(Mcmc.Post.mean, Symbol(name))[2:end,:])
    else
        True = getfield(Mcmc.truePara, Symbol(name))
        Esti = getfield(Mcmc.Post.mean, Symbol(name))
    end
    Diff = abs.(Esti .- True)
    
    names = ["Esti", "True", "|Diff|"]
    values = [Esti True Diff]
    values = round.(values, digits=3)


    # print title
    println(join(names, "\t"))
    println(join(fill("=======", length(names)), "\t"))
    # print data
    for row in 1:size(values, 1)     
        rounded_row = round.(values[row,:], digits=3)
        println(join(rounded_row, "\t"))
    end
end


# ╔═╡ 8a35299b-8fd9-41d3-8ed1-a88312d1ba43
"""
"""
function checkConvergence(MCMC)

    mcmcRa = MCMC.Post.ra[(MCMC.Cond.nBurnin+1):end, :, :]
    essRa = Chains(mcmcRa) |> ess_rhat

    mcmcRt = MCMC.Post.rt[(MCMC.Cond.nBurnin+1):end, :, :]
    essRt = Chains(mcmcRt) |> ess_rhat

    mcmcQr = MCMC.Post.qr[(MCMC.Cond.nBurnin+1):end, :, :]
    essQr = Chains(mcmcQr) |> ess_rhat

    essLength = length(essRa[:,2]) + length(essRt[:,2]) + count(!isnan, (essQr[:,2]))
    rhatLength = length(essRa[:,3]) + length(essRt[:,3]) + count(!isnan, vec(essQr[:,3]))
    essOk = count(essRa[:,2] .> 400) + count(essRt[:,2] .> 400) + count( skipmissing(essQr[:,2] .> 400) )
    rhatOk = count(essRa[:,3] .< 1.1) + count(essRt[:,3] .< 1.1) + count( skipmissing(essQr[:,3] .< 1.1) )

    Conv = (
        ess = (essOk / essLength)*100,
        rhat = (rhatOk / rhatLength)*100,
        essN = "$(essOk) / $(essLength)",
        rhatN = "$(rhatOk) / $(rhatLength)"
        
    )
    return Conv
end

# ╔═╡ 95258e82-f774-43c3-8056-62c935512a7c

## Start a one-condition Simulation Study!
"""
    runSimulation(Cond, truePara; Para=(:a, :b, :λ, :σ²t), funcData=setDataRtIrt, funcGibbs=GibbsRtIrt)

### Arguments.
	- Para. can be any name in the `InputPara`-struct (:a, :b, :λ, :σ²t, :ρ, :β, :Σp).
	- funcData.
	- funcGibbs.

"""
function runSimulation(Cond, truePara; Para=(:a, :b, :λ, :σ²t), funcData=setDataRtIrt, funcGibbs=GibbsRtIrt, typeName="norm")
    fooData = getfield(Main, Symbol(funcData))
    fooGibbs = getfield(Main, Symbol(funcGibbs))

    Run = Dict(Symbol(:True) => [], [run => Dict() for run in 1:Cond.nRep]...)
    True = Dict(p => [] for p in Para)
    for p in Para
        True[p] = vec(getfield(truePara, p))
    end
    Run[:True] = True

    for run in 1:Cond.nRep

        # 在中間位置加入清理
        if mod(run, Cond.nRep/10) == 0
            GC.gc(true)
        end

        ## Data.
        Data = fooData(Cond, truePara; type=typeName)

        ## Fit.
        MCMC = fooGibbs(Cond, truePara=truePara, Data=Data)
        sample!(MCMC)   

        ## Save Post Means and Dic.
        Post = Dict(Symbol(:Dic) => [], Symbol(:Diag) => [], [p => [] for p in Para]...)
        for p in Para
            Post[p] = getfield(MCMC.Post.mean, p)
        end
        Post[:Dic] = [getDic(MCMC).DIC]
		Post[:Diag] = checkConvergence(MCMC)
        Run[run] = Post
    end
    return Run
end

# ╔═╡ 14a223bc-4419-4773-94b8-82aef171c285
md"""

### A Basic Example about Conducting a Simulation.

```

using ExtendedRtIrtModeling,
    Random,
    Plots,
    TidierPlots,
    DataFrames,
    JLD2

## ==========================
##  Just for testing.
## ==========================


## 30104
begin
## Conditions.
Random.seed!(1234)

## Manipulate the conditions
# df30104 = Dict("$(J)-$(K)" => Dict() for J in (250, 500, 1000), K in (15) )
for J in (250, 500, 1000), K in (15)
    @time begin 

        # 1. 設定條件
        Cond = setCond(nSubj=J, nItem=K, nRep=100, nIter=10_000, nChain=3)
        truePara = setTrueParaRtIrt(Cond)
        
        # 2. 執行模擬
        res = runSimulation(Cond, truePara; 
                           funcData=setDataRtIrt, 
                           funcGibbs=GibbsRtIrt, 
                           Para=(:a, :b, :λ, :σ²t, :β, :Σp))

        # 3. 立即儲存這個條件的結果
        condName = "$(J)-$(K)"
        @save "sim30104_$(condName).jld2" res

        # 4. 清理記憶體
        GC.gc()

        println("===== (=^・^=) ===== Condition $(condName) DONE! ===== (=^・^=) =====")
    
    end    
end


# 5. 最後如果需要合併所有結果
df30104 = Dict()
for J in (250, 500, 1000), K in (15)
    condName = "$(J)-$(K)"
    @load "sim30104_$(condName).jld2" res
    df30104[condName] = res
end

# 6. 儲存完整結果
@save "sim30104_full.jld2" df30104

end


```


"""

# ╔═╡ b8a44521-8622-4b83-9518-a2da3107d91e
"""
function evaluate(MCMC::GibbsRtIrtNull)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end
"""

# ╔═╡ cf9fd679-b387-476f-8401-46082d4f6c93
"""
function evaluate(MCMC::GibbsRtIrtNull)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end
"""

# ╔═╡ 8d47ab56-b591-4d18-85c2-03550463d90a
"""
function evaluate(MCMC::GibbsRtIrtQuantile)
	Rmse = [getRmse(MCMC.Post.mean.a , MCMC.truePara.a) getRmse(MCMC.Post.mean.b , MCMC.truePara.b) getRmse(MCMC.Post.mean.β, MCMC.truePara.β)]
	Bias = [getBias(MCMC.Post.mean.a, MCMC.truePara.a) getBias(MCMC.Post.mean.b , MCMC.truePara.b) getBias(MCMC.Post.mean.β, MCMC.truePara.β)]
	Dic = getDicMlIrt(MCMC.Cond, MCMC.Data, MCMC.Post).DIC
	return Rmse, Bias, Dic
end
"""

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
# ╟─3dd409f0-ab34-11ef-22d1-3b72e5416118
# ╠═1d4a9f51-3db8-4a78-936d-a7f2879cc0a0
# ╠═23484455-73cd-4c7f-a13d-f58fb0bc3051
# ╠═bc0d18af-ab68-443d-8e98-ec5fbedc2e12
# ╠═b57756d7-b28a-4fba-be7d-d1b9f196fb20
# ╠═d66e453f-de50-42a5-a821-508e23b44b10
# ╠═3468403a-3ae8-4635-9fd9-66538be96315
# ╠═74c8ac72-0701-4a5a-8409-0d93e6022503
# ╠═201f594a-cc89-4236-957f-ed88953b9ebd
# ╠═595dbac9-1bb6-4975-928c-5a53f28eae56
# ╠═838de637-d131-480b-99d8-ed55433f41de
# ╠═e3e68b65-9e25-49c7-8a8f-07b34dfa1b0f
# ╠═732cc1f5-549e-492f-a37d-86ed70cdb7b8
# ╠═f8e59e9a-16d6-4e4e-b630-f81f593d9f01
# ╠═54a39e13-807d-49ba-9cd0-710c7e500eff
# ╠═7fa851f0-f16d-43d1-af39-a992f5764c88
# ╠═7b1041d1-bb69-412a-80e8-c0ab731545f4
# ╠═cf125fa8-18c0-4f44-80fd-4447c45a5876
# ╠═d792cfdb-503f-4986-a62d-02338db36ac4
# ╠═95258e82-f774-43c3-8056-62c935512a7c
# ╠═8a35299b-8fd9-41d3-8ed1-a88312d1ba43
# ╠═14a223bc-4419-4773-94b8-82aef171c285
# ╠═b8a44521-8622-4b83-9518-a2da3107d91e
# ╠═cf9fd679-b387-476f-8401-46082d4f6c93
# ╠═8d47ab56-b591-4d18-85c2-03550463d90a
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

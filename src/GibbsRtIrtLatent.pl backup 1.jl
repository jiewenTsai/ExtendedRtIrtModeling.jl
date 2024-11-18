### A Pluto.jl notebook ###
# v0.19.46

using Markdown
using InteractiveUtils

# ╔═╡ a9b62015-e7d7-452e-8673-4c5984274422
using PlutoUI,
	ProgressMeter

# ╔═╡ 6ef8b4d0-9f55-11ef-3908-19b16335c1cc
md"""
# GibbsRtIrtLatent

Including 2 classes,

- GibbsRtIrtLatent
- GibbsRtIrtLatentQr


"""

# ╔═╡ 87798f0b-fb5f-43ef-9205-054961dcdc6b
TableOfContents()

# ╔═╡ b0c58544-da9c-4aaa-b96e-d21b3c22e7df
md"""
## Structs
"""

# ╔═╡ acc558d5-3644-4b5a-b9ca-09e62260c194
mutable struct OutputPostRtIrtLatent
	ra::Array
	rt::Array
	qr::Array
	logLike::Array
	mean
	function OutputPostRtIrtLatent(Cond; ra=[], rt=[], qr=[], logLike=[], mean=Float64[]
)
		ra=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain)
		rt=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain,)
		qr=Array{Float64}(undef, Cond.nIter, (Cond.nFeat+2)+4, Cond.nChain)
		logLike=Array{Float64}(undef, Cond.nIter, 1, Cond.nChain)
		return new(ra, rt, qr, logLike, mean)
	end
end

# ╔═╡ 8da2f52b-b034-4b24-9216-f19b9729db70
mutable struct OutputPostRtIrtLatentQr
	ra::Array
	rt::Array
	qr::Array
	logLike::Array
	mean
	function OutputPostRtIrtLatentQr(Cond; ra=[], rt=[], qr=[], logLike=[], mean=Float64[]
)
		ra=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain)
		rt=Array{Float64}(undef, Cond.nIter, Cond.nSubj+Cond.nItem+Cond.nItem, Cond.nChain,)
		qr=Array{Float64}(undef, Cond.nIter, (Cond.nFeat+2+4+Cond.nSubj), Cond.nChain)
		logLike=Array{Float64}(undef, Cond.nIter, 1, Cond.nChain)
		return new(ra, rt, qr, logLike, mean)
	end
end

# ╔═╡ 8b51811b-490f-47fb-91c7-84c83eff3dca
mutable struct GibbsRtIrtLatent
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsRtIrtLatent)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			ζ=randn(self.Cond.nSubj),
			λ = zeros(self.Cond.nItem),
			σ²t = ones(self.Cond.nItem),
			β=randn(self.Cond.nFeat+2),
            Σp=I(2)  )
		return self
	end
	function GibbsRtIrtLatent(Cond; 
	Data = Float64[], 
	truePara = Float64[], 
	Para = Float64[],
	Post = Float64[]
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPostRtIrtLatent(Cond)
		return obj
	end

end

# ╔═╡ d531a3ed-d0d0-45c5-9094-abc6bfadf31e
mutable struct GibbsRtIrtLatentQr
	Cond
	Data
	truePara
	Para
	Post
	## Initiation

	function setInitialValues(self::GibbsRtIrtLatentQr)
		self.Para = InputPara(
			θ=randn(self.Cond.nSubj), 
			a=ones(self.Cond.nItem), 
			b=zeros(self.Cond.nItem), 
			ζ=randn(self.Cond.nSubj),
			λ = zeros(self.Cond.nItem),
			σ²t = ones(self.Cond.nItem),
			β=randn(self.Cond.nFeat+2),
            Σp=I(2)  )
		return self
	end
	function GibbsRtIrtLatentQr(Cond; 
	Data = Float64[], 
	truePara = Float64[], 
	Para = Float64[],
	Post = Float64[]
	)			
		obj = new(Cond, Data, truePara, Para, Post)
		setInitialValues(obj)
		obj.Post = OutputPostRtIrtLatentQr(Cond)
		return obj
	end

end

# ╔═╡ ad214791-4576-440b-9c21-99d4c41b4a01
md"""
## sample

"""

# ╔═╡ 25c6b329-28c2-460f-b6db-b532d8918f64

"""
	sample!(MCMC::GibbsRtIrtLatent)
"""
function sample!(MCMC::GibbsRtIrtLatent; intercept=false, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	
		## structural
		#Para.ν = drawQrWeightsLatentQr(Cond,Data,Para)
		Para.β = getSubjCoefficientsLatent(Cond,Data,Para)
		#Para.β = zeros(Cond.nFeat+1,2)
		if intercept==false
			Para.β[1] = 0.
		end
		Para.Σp = drawSubjCovarianceLatent(Cond, Data, Para)
		
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
		Para.ζ = drawSubjSpeedLatent(Cond,Data,Para)

		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.ζ; Para.λ; Para.σ²t]
		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodRtIrtLatent(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(

		## structure
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:((Cond.nFeat+2)), :], dims=(1,3)) |> vec,
		Σp = mean(Post.qr[(Cond.nBurnin+1):end, ((Cond.nFeat+2)+1):end, :], dims=(1,3)) |> vec,

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

# ╔═╡ 070dd7b7-7463-45db-b06e-a22caee0935a

"""
	sample!(MCMC::GibbsRtIrtLatentQr)
"""
function sample!(MCMC::GibbsRtIrtLatentQr; intercept=false, itemtype::Union{String} = "2pl")

	if !(itemtype in ["1pl", "2pl"])
        error("Invalid input: the item type must be '1pl' or '2pl'.")
    end
	
	Cond = MCMC.Cond
	Data = MCMC.Data
	Para = MCMC.Para
	Post = MCMC.Post

	@showprogress for m in 1:Cond.nIter, l in 1:Cond.nChain	
		## structural
		Para.ν = drawQrWeightsLatentQr(Cond,Data,Para)
		Para.β = getSubjCoefficientsLatentQr(Cond,Data,Para)
		#Para.β = zeros(Cond.nFeat+1,2)
		if intercept==false
			Para.β[1] = 0.
		end
		Para.Σp = I(2)  #drawSubjCovarianceLatentQr(Cond, Data, Para)
		
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
		Para.ζ = drawSubjSpeedLatentQr(Cond,Data,Para)

		## collect
		Post.qr[m, :, l] = [vec(Para.β); vec(Para.Σp); vec(Para.ν)]
		Post.ra[m, :, l] = [Para.θ; Para.a; Para.b]
		Post.rt[m, :, l] = [Para.ζ; Para.λ; Para.σ²t]
		## logLike
		Post.logLike[m, :, l] .= getLogLikelihoodRtIrtLatentQr(Cond,Data,P=Para)
	end

	## summarystats
	Post.mean = InputPara(

		## structure
		β = mean(Post.qr[(Cond.nBurnin+1):end, 1:((Cond.nFeat+2)), :], dims=(1,3)) |> vec,
		Σp = mean(Post.qr[(Cond.nBurnin+1):end, (Cond.nFeat+2+1):(Cond.nFeat+2+4), :], dims=(1,3)) |> vec,
		ν = mean(Post.qr[(Cond.nBurnin+1):end, (Cond.nFeat+2+4+1):end, :], dims=(1,3)) |> vec,

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

# ╔═╡ ab3a0f30-ce07-4045-873a-c238941d4092
md"""
## Functions
"""

# ╔═╡ 09ed579c-de74-4010-82e5-374da92cc565
function getLogLikelihoodRtIrtLatent(Cond,Data; P=Para)
    x = [ones(Cond.nSubj) Data.X P.θ]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.λ' .- P.ζ
    μOfζ = x * P.β
	Σp = reshape(P.Σp,2,2)
    
    logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, sqrt.(P.σ²t')), Data.logT)) sum(logpdf.(Normal.(μOfζ, sqrt.(Σp[2,2])), P.ζ)) ] 

    return sum(logPdf)
end

# ╔═╡ 1a5da9da-2660-4f57-b483-2ec6b856fc1a
md"""
!!! info "Think about the Σp"
	
	
"""

# ╔═╡ 8ea32a2a-f43b-482b-837a-bd2b94a903c3
function getLogLikelihoodRtIrtLatentQr(Cond,Data; P=Para)

	k1Rt = (1 - 2 * Cond.qRt) / (Cond.qRt * (1 - Cond.qRt))
    k2Rt = 2 / (Cond.qRt * (1 - Cond.qRt))
	eRt = reshape(P.ν, Cond.nSubj)
    k1e = k1Rt.*eRt
	k2e = (k2Rt.*eRt)

	
    x = [ones(Cond.nSubj) Data.X P.θ]
    pr = P.a' .* (P.θ .- P.b')
    μt = P.λ' .- P.ζ
    
	#Σp = reshape(P.Σp,2,2)
	μOfζ = x * P.β .+ k1e 

	#sum([logpdf.(Normal.(μOfζ[i,:], sqrt.(Σp[2,2] .* k2e[i,:])), P.ζ[i,:]) for i in 1:Cond.nSubj])
    logPdf = [sum(logpdf.(BernoulliLogit.(pr), Data.Y)) sum(logpdf.(Normal.(μt, sqrt.(P.σ²t')), Data.logT)) sum(logpdf.(Normal.(μOfζ, sqrt.(1)), P.ζ)) ] 

    return sum(logPdf)
end

# ╔═╡ 0cff6c79-e493-46b8-a638-113b6e675869
"""
"""
function getDic(MCMC::GibbsRtIrtLatent)
	DIC = OutputDic()	
	## compute dic
	D̂ = -2*getLogLikelihoodRtIrtLatent(MCMC.Cond, MCMC.Data, P= MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD
	return DIC
	
end

# ╔═╡ d96a81d5-c364-4314-a3cf-be66f548e566
"""
"""
function getDic(MCMC::GibbsRtIrtLatentQr)
	DIC = OutputDic()	
	## compute dic
	D̂ = -2*getLogLikelihoodRtIrtLatentQr(MCMC.Cond, MCMC.Data, P= MCMC.Post.mean)
	D̄ = -2*mean(MCMC.Post.logLike)
	DIC.pD = D̄ - D̂
	DIC.DIC = D̄ + DIC.pD
	return DIC
	
end

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
ProgressMeter = "92933f4c-e287-5a05-a399-4b506db050ca"

[compat]
PlutoUI = "~0.7.60"
ProgressMeter = "~1.10.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.11.0"
manifest_format = "2.0"
project_hash = "5e7f29b122549dc0b0bee5281e0d2ec115e16bd5"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "6e1d2a35f2f90a4bc7c2ed98079b2ba09c35b83a"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.2"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.2"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.11.0"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"
version = "1.11.0"

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
version = "1.11.0"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"
version = "1.11.0"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"
version = "1.11.0"

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
version = "1.11.0"

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
version = "8.6.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"
version = "1.11.0"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.7.2+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"
version = "1.11.0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
version = "1.11.0"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"
version = "1.11.0"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"
version = "1.11.0"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.6+0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"
version = "1.11.0"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.12.12"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.27+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "Random", "SHA", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.11.0"

    [deps.Pkg.extensions]
    REPLExt = "REPL"

    [deps.Pkg.weakdeps]
    REPL = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

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
version = "1.11.0"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "8f6bc219586aef8baf0ff9a5fe16ee9c70cb65e4"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.10.2"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"
version = "1.11.0"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"
version = "1.11.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"
version = "1.11.0"

[[deps.Statistics]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "ae3bb1eb3bba077cd276bc5cfc337cc65c3075c0"
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.11.1"

    [deps.Statistics.extensions]
    SparseArraysExt = ["SparseArrays"]

    [deps.Statistics.weakdeps]
    SparseArrays = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

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
version = "1.11.0"

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
version = "1.11.0"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"
version = "1.11.0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.11.0+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.59.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"
"""

# ╔═╡ Cell order:
# ╟─6ef8b4d0-9f55-11ef-3908-19b16335c1cc
# ╠═a9b62015-e7d7-452e-8673-4c5984274422
# ╠═87798f0b-fb5f-43ef-9205-054961dcdc6b
# ╠═b0c58544-da9c-4aaa-b96e-d21b3c22e7df
# ╠═acc558d5-3644-4b5a-b9ca-09e62260c194
# ╠═8da2f52b-b034-4b24-9216-f19b9729db70
# ╠═8b51811b-490f-47fb-91c7-84c83eff3dca
# ╠═d531a3ed-d0d0-45c5-9094-abc6bfadf31e
# ╟─ad214791-4576-440b-9c21-99d4c41b4a01
# ╠═25c6b329-28c2-460f-b6db-b532d8918f64
# ╠═070dd7b7-7463-45db-b06e-a22caee0935a
# ╠═ab3a0f30-ce07-4045-873a-c238941d4092
# ╠═09ed579c-de74-4010-82e5-374da92cc565
# ╟─1a5da9da-2660-4f57-b483-2ec6b856fc1a
# ╠═8ea32a2a-f43b-482b-837a-bd2b94a903c3
# ╠═0cff6c79-e493-46b8-a638-113b6e675869
# ╠═d96a81d5-c364-4314-a3cf-be66f548e566
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

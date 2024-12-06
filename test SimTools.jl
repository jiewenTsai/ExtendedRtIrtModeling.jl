
include("src/ExtendedRtIrtModeling.jl")
import .ExtendedRtIrtModeling as ex

using Distributions,
    LinearAlgebra,
    Plots,
    StatsPlots,
    HTTP,
    CSV,
    DataFrames,
    Random,
    JLD2



Random.seed!(1234)
Cond = ex.setCond(nChain=2, nIter=3000, nSubj=500, nItem=15, qRt=0.25)
True = ex.setTrueParaRtIrtCross(Cond)
#Data = ex.setDataRtIrt(Cond, True)
Data = ex.setDataRtIrtCross(Cond, True, type="skew")
Mcmc = ex.GibbsRtIrtCross(Cond, Data=Data, truePara=True)
ex.sample!(Mcmc)


DataNorm = ex.setDataRtIrtCross(Cond, True, type="norm")
DataSkew = ex.setDataRtIrtCross(Cond, True, type="skew")
DataTail = ex.setDataRtIrtCross(Cond, True, type="tail")

density( vec(DataNorm.logT) )
density( vec(DataSkew.logT))
density( (DataTail.logT))



rand(LogNormal(0, sqrt.(True.σ²t')), Cond.nSubj)

ex.coef(Mcmc)
ex.comparePara(Mcmc, name=:σ²t )
ex.comparePara(Mcmc, name=:λ )
ex.comparePara(Mcmc, name=:ρ )
ex.comparePara(Mcmc, name=:b )
ex.comparePara(Mcmc, name=:β )

Mcmc.Para.ν

comboJK02 = [(1000,30)]
comboQT25 = [(0.25, "norm"), (0.25, "skew")]

for (J,K) in comboJK02, (Q,Type) in comboQT25
    #try
        @time begin 
            # 1. 設定條件
            Cond = ex.setCond(nSubj=J, nItem=K, nRep=10, nIter=2_000, nChain=2, qRt=Q)
            True = ex.setTrueParaRtIrtCross(Cond)
            
            # 2. 執行模擬
            Data = ex.setDataRtIrtCross(Cond, True, type=Type)
            Mcmc = ex.GibbsRtIrtCrossQr(Cond, Data=Data, truePara=True)
            ex.sample!(Mcmc)

            """
            res = ex.runSimulation(Cond, truePara; 
                            funcData=setDataRtIrtCross, 
                            funcGibbs=GibbsRtIrtCrossQr, 
                            typeName=Type,
                            Para=(:a, :b, :λ, :σ²t, :ρ, :Σp))
            """
            # 3. 立即儲存這個條件的結果
            condName = "$(J)-$(K)-$(Q)-$(Type)"
            #@save "sim30200_$(condName).jld2" res

            # 4. 清理記憶體
            #res = nothing
            GC.gc(true)

            println("===== (=^.^=) ===== Condition $(condName) DONE! ===== (=^.^=) =====")
        
        end    
    """ 
    catch e
        condName = "$(J)-$(K)-$(Q)-$(Type)"
        println("跳過條件 $(condName) 發生錯誤: ", e)
        continue  # 直接進行下一次迭代
    end
    """
end





ex.coef(Mcmc)





























Cond = ex.setCond(nChain=2, nIter=3000, nSubj=500, nItem=10)
truePara = ex.setTrueParaRtIrt(Cond, trueCorr=0.5)
Data = ex.setDataRtIrtNull(Cond, truePara)
#Data2 = ex.setDataMlIrt(Cond, truePara)

## generate data sets

COND = ex.setCond(nChain=2, nIter=3000, nSubj=1000, nItem=10, nRep=10)
truePARA = ex.setTrueParaRtIrt(COND)

#begin
data_list = Array{Float64}(undef, COND.nSubj, COND.nItem+COND.nItem+COND.nFeat, COND.nRep)
for i in 1:COND.nRep
    DATA = ex.setDataRtIrt(COND, truePARA)
    data_list[:,:,i] = [DATA.Y DATA.T DATA.X]
end

data_list

using JSON

json_data = JSON.json(data_list)
write("datasets10.json", json_data)
#end

MCMC = ex.GibbsRtIrtCross(Cond, Data=Data, truePara=truePara)
ex.sample!(MCMC)

ex.coef(MCMC)


ex.getRmseBasic(MCMC.truePara.a, MCMC.Post.mean.a)
ex.getBias(MCMC.truePara.a, MCMC.Post.mean.a)
cor(MCMC.truePara.a, MCMC.Post.mean.a)
[MCMC.truePara.a MCMC.Post.mean.a]
ex.getRmse(MCMC.truePara.b, MCMC.Post.mean.b)
ex.getBias(MCMC.truePara.b, MCMC.Post.mean.b)
cor(MCMC.truePara.b, MCMC.Post.mean.b)
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

sqrt(sum((MCMC.Post.mean.a .- MCMC.truePara.a).^2) / 10)
sqrt(mean((MCMC.Post.mean.b .- MCMC.truePara.b).^2))

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



D = ex.testingDict(300,10,5)
D["T"]


# =======
# sim data
# ========


Random.seed!(1234)
Cond = ex.setCond(nSubj=10, nItem=3)
truePara = ex.setTrueParaRtIrt(Cond)
Data = ex.setDataRtIrt(Cond, truePara)
Data.T

## generate data sets
Random.seed!(1234)
COND = ex.setCond(nSubj=10, nItem=3)
truePARA = ex.setTrueParaRtIrt(Cond)

data_list = Array{Float64}(undef, COND.nSubj, COND.nItem+COND.nFeat, COND.nRep)
for i in 1:COND.nRep
    DATA = ex.setDataRtIrt(COND, truePARA)
    data_list[:,:,i] = [DATA.Y DATA.X]
end

data_list




@progress for n in 1:COND.nRep
 
    DATA = setData(COND, truePARA)
    MCMC = sample!(GibbsSamplerPgMlIrt(DATA=DATA, truePARA=truePARA))
    EVAL.Rmse[:,n],EVAL.Bias[:,n],EVAL.Dic[n] = evaluate(MCMC)
    

end


str = "MCMC.Post.mean." .* "a"

eval(Symbol(str))

getRmse2(name) = mean(sqrt(mean((Base.getproperty(MCMC.Post.mean, name) .- Base.getproperty(MCMC.truePara, name)).^2)))

# ======================

struct OutputMetrics
    Rmse 
    Bias 
    Corr 
    Dic 
    function OutputMetrics(
        Rmse = DataFrame(),
        Bias = DataFrame(),
        Corr = DataFrame(),
        Dic = []
    )
        return new(Rmse, Bias, Corr, Dic)
    end
end

## ========================================================
## running a simulation study.
## ========================================================

## Conditions.
Random.seed!(1234)
Cond = ex.setCond(nSubj=100, nItem=3, nRep=3)
truePara = ex.setTrueParaRtIrt(Cond)

## data container 
df = OutputMetrics()

## Start Simulation Study!
for run in 1:Cond.nRep
    dictRmse = Dict()
    dictBias = Dict()
    dictCorr = Dict()
    arrayDic = Float64[]

    ## Data.
    Data = ex.setDataRtIrt(Cond, truePara)

    ## Fit.
    Mcmc = ex.GibbsRtIrt(Cond, truePara=truePara, Data=Data)
    ex.sample!(Mcmc)

    ## Save Metrics.
    for i in (:a, :b, :λ, :σ²t)
        metricRmse = ex.getRmse(Mcmc, i) 
        metricBias = ex.getBias(Mcmc, i) 
        metricCorr = ex.getCorr(Mcmc, i) 
        dictRmse[i] = metricRmse
        dictBias[i] = metricBias
        dictCorr[i] = metricCorr
    end
    push!(df.Rmse, (run=run, dictRmse...)) 
    push!(df.Bias, (run=run, dictBias...))
    push!(df.Corr, (run=run, dictCorr...))

    ### Dic
    metricDic = ex.getDic(Mcmc).DIC
    push!(df.Dic, metricDic)
end

df.Dic

@save "sim3.jld2" df

@load "sim3.jld2"

df2 = load("sim3.jld2")

df2["df"].Rmse

# 轉換為 DataFrame
df = DataFrame()
for (run, mm) in res
    row = merge(Dict(:run => run), mm)  # 加入運算次數
    push!(df, row)
end


ex.getBias(MCMC, :a) 
ex.getCorr(MCMC, :a) 

Base.getproperty(MCMC.Post.mean, :a)


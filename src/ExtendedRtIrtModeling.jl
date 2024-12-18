module ExtendedRtIrtModeling

#include("GenInvGaussian.jl")
#import .GenInvGaussian as ge

using Markdown
    #InteractiveUtils

using ProgressMeter,
	ProgressLogging,
    #SimplePlutoInclude,
    DataFrames,
    LinearAlgebra,
    Kronecker,
    PolyaGammaSamplers,
    Distributions,
    Random,
    Plots,
    MCMCChains,
    PrettyTables,
    PlutoUI


#include("notebook1102.jl")
include("Base.pl.jl")
include("Draw.pl.jl")
include("SimTools.jl")
#include("GibbsMlIrt.pl.jl")
include("GibbsRtIrt.pl.jl")
include("GibbsRtIrtCross.pl.jl")
include("GibbsRtIrtLatent.pl.jl")

export 
	## set
	setCond, 
    setData, 
    InputData, 
    setDataMlIrt, 
    setDataRtIrt, 
    setDataRtIrtNull,
    setDataRtIrtCross,
    setDataRtIrtLatent,
	setTrueParaMlIrt, 
    setTrueParaRtIrt,
    setTrueParaRtIrtCross,
    setTrueParaRtIrtLatent,

	## get
	getBias, 
    getRmse, 
    getCorr,
    getDic,
    getMetrics,
    getMetrics2,
    getRound3,
	#evaluate,

	## sample!
	sample!,

	## structs
	GibbsMlIrt, 
    GibbsRtIrtNull, 
    GibbsRtIrt, 
    #GibbsRtIrtQuantile,
    GibbsRtIrtCross,
    GibbsRtIrtCrossQr,
    GibbsRtIrtLatent,
    GibbsRtIrtLatentQr,
	
	## Miscellaneous
	coef, 
    precis,
    comparePara,
    checkConvergence,
    OutputMetrics,
    runSimulation    



end
module PearsonDistribution
using Distributions
using QuadGK
using SpecialFunctions


import Base.convert
import Distributions: params, partype, pdf, logpdf, logcdf, quantile, rand, cdf, 
       @distr_support, @check_args, minimum, maximum

include("pearson.jl")
include("pearsonIV.jl")
include("fitpearson.jl")

export
PearsonI, PearsonII, PearsonIII, PearsonIV, PearsonV, PearsonVI, PearsonVII,
fitpearson
end # module

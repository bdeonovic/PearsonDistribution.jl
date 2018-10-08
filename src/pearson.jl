const PearsonI{T} = LocationScale{T, Distributions.Beta{T}} where T <: Real
PearsonI(a::T, b::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, Beta{T}}(location, scale, Beta(a, b))
PearsonI(a::Real, b::Real, location::Real, scale::Real) = 
  PearsonI(promote(a, b, location, scale)...)
PearsonI(a::Integer, b::Integer, location::Integer, scale::Integer) = 
  PearsonI(Float64(a), Float64(b), Float64(location), Float64(scale))

const PearsonII{T} = LocationScale{T, Distributions.Beta{T}} where T <: Real
PearsonII(a::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, Beta{T}}(location, scale, Beta(a, a))
PearsonII(a::Real, location::Real, scale::Real) = 
  PearsonII(promote(a, location, scale)...)
PearsonII(a::Integer, location::Integer, scale::Integer) = 
  PearsonII(Float64(a), Float64(location), Float64(scale))

const PearsonIII{T} = LocationScale{T, Distributions.Gamma{T}} where T <: Real
PearsonIII(shape::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, Gamma{T}}(location, one(T), Gamma(shape, scale))
PearsonIII(shape::Real, location::Real, scale::Real) = 
  PearsonIII(promote(shape, location, scale)...)
PearsonIII(shape::Integer, location::Integer, scale::Integer) = 
  PearsonIII(Float64(shape), Float64(location), Float64(scale))

const PearsonV{T} = LocationScale{T, Distributions.InverseGamma{T}} where T <: Real
PearsonV(shape::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, InverseGamma{T}}(location, one(T), InverseGamma(shape, scale))
PearsonV(shape::Real, location::Real, scale::Real) = 
  PearsonV(promote(shape, location, scale)...)
PearsonV(shape::Integer, location::Integer, scale::Integer) = 
  PearsonV(Float64(shape), Float64(location), Float64(scale))

const PearsonVI{T} = LocationScale{T, Distributions.BetaPrime{T}} where T <: Real
PearsonVI(a::T, b::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, BetaPrime{T}}(location, scale, BetaPrime(a, b))
PearsonVI(a::Real, b::Real, location::Real, scale::Real) = 
  PearsonVI(promote(a, b, location, scale)...)
PearsonVI(a::Integer, b::Integer, location::Integer, scale::Integer) = 
  PearsonVI(Float64(a), Float64(b), Float64(location), Float64(scale))

const PearsonVII{T} = LocationScale{T, Distributions.TDist{T}} where T <: Real
PearsonVII(df::T, location::T, scale::T) where T <: Real = 
  LocationScale{T, TDist{T}}(location, scale, TDist(df))
PearsonVII(df::Real, location::Real, scale::Real) = 
  PearsonVII(promote(df, location, scale)...)
PearsonVII(df::Integer, location::Integer, scale::Integer) = 
  PearsonVII(Float64(df), Float64(location), Float64(scale))

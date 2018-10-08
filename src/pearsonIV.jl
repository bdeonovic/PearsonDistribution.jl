struct PearsonIV{T<:Real} <: ContinuousUnivariateDistribution
  m::T
  nu::T
  location::T
  scale::T

  function PearsonIV{T}(m::T, nu::T, location::T, scale::T) where T <: Real
    @check_args(PearsonIV, m > 0.5 && scale > zero(0.0))
    new{T}(m, nu, location, scale)
  end
end

PearsonIV(m::T, nu::T, location::T, scale::T) where T <: Real = 
  PearsonIV{T}(m, nu, location, scale)
PearsonIV(m::Real, nu::Real, location::Real, scale::Real) = 
  PearsonIV(promote(m, nu, location, scale)...)
function PearsonIV(m::Integer, nu::Integer, location::Integer, scale::Integer)
  PearsonIV(Float64(m), Float64(nu), Float64(location), Float64(scale))
end
PearsonIV(m::Real, nu::Real) = PearsonIV(m, nu, 0.0, 1.0)

@distr_support PearsonIV -Inf Inf

#### Conversions
function convert(::Type{PearsonIV{T}}, m::Real, nu::Real, location::Real, 
                 scale::Real) where T <: Real
  PearsonIV(T(m), T(nu), T(location), T(scale))
end
function convert(::Type{PearsonIV{T}}, d::PearsonIV{S}) where{T <: Real, S <: Real}
  PearsonIV(T(d.m), T(d.nu), T(d.location), T(d.scale))
end

#### Parameters

params(d::PearsonIV) = (d.m, d.nu, d.location, d.scale)
@inline partype(d::PearsonIV{T}) where {T<:Real} = T

#### Functions
function gammar2(x::T, y::T) where T <: Real
  #returns abs(gamma(x+iy)/gamma(x))^2
  y2 = y^2
  xmin = (2*y2 > 10.0) ? 2*y2 : 10.0
  r=1
  s=1
  p=1
  f=0
  while x < xmin
    t = y/x
    x+=1
    r *= (1 + t^2)
  end
  while p > (s * eps())
    p *= y2 + f^2
    f += 1
    p /= x * f
    x += 1
    s += p
  end
  return 1.0 / (r * s)
end

function logP4norm(m::T, nu::T, a::T) where T <: Real
  logk = - 0.5 * log(pi) - log(a) - lgamma(m - 0.5) + 
    lgamma(m) + log(gammar2(m, 0.5 * nu))
  return logk
end

function logpdf(d::PearsonIV, x::Real)
  logk = logP4norm(d.m, d.nu, d.scale)
  return logk - d.m * log(1 + ((x - d.location) / d.scale) ^ 2) - 
    d.nu * atan((x - d.location) / d.scale)
end

pdf(d::PearsonIV, x::Real) = exp(logpdf(d, x))

function rand(d::PearsonIV)
  if d.m < 1 
    error("Simulation only for m >= 1")
  end
  logk = logP4norm(d.m, d.nu, d.scale)
  return rand_logpearsonIV_logk(d.m, d.nu, d.location, d.scale, logk)
end

function rand_logpearsonIV_logk(m::T, nu::T, location::T, scale::T, logk::T) where T <: Real
  b = 2 * m - 2
  Mode = atan(-nu/b)
  cosM = b / sqrt(b ^ 2 + nu ^ 2)
  r = b * log(cosM) - nu * Mode
  rc = exp(-r - logk)

  x = 0.0
  z = 0.0
  s = 0

  iter = 0

  while (iter == 0) || (abs(x) >= (pi / 2)) || 
        ((z + log(rand())) > (b * log(cos(x)) - nu * x - r))
    s = 0
    z = 0
    x = 4 * rand()
    if x > 2
      x -= 2
      s = 1
    end
    if x > 1
      z = log(x - 1)
      x = 1 - z
    end
    x = s != 0 ? Mode + rc * x : Mode - rc * x
    iter += 1
  end
  return scale * tan(x) + location
end

function cdf(d::PearsonIV, x::Real, tol::Float64 = 1.0e-8)
  return cdf_pearson4_num_int(x, d.m, d.nu, d.location, d.scale, tol)
end

function cdf_pearson4_num_int(x::Real, m::T, nu::T, location::T, scale::T, tol::Float64) where {T<:Real}
  modus = location - scale * nu / (2 * m)
  logk = logP4norm(m, nu, scale)
  f(z) = exp(logk - m * log(1 + ((z - location) / scale)^2) - nu * atan((z - location) / scale))
  res = 0.0
  tol = 0.0
  if x > 0
    res, atol = quadgk(f, x, Inf, rtol=tol)
    res = 1 - res
  else
    res, atol = quadgk(f, -Inf, x, rtol=tol)
  end
  return res
end


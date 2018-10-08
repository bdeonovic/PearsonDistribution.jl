function fitpearson(mu1::T, mu2::T, beta1::T, beta2::T) where T <: Real
  if beta1 > (beta2 - 1)
    error("There are no probability distributions with these moments")
  elseif isapprox(beta1, (beta2 -1))
    p = (1 + sqrt(beta1)/sqrt(4+beta1))/2
    a = mu1 - mu2 * sqrt((1-p)/p)
    b = mu1 + mu2 * sqrt(p/(1-p))
    error("Not a pearson distribution. Try discrete distribution with mass \n $p on $a and mass $(1-p) on $b")
  else
    b0 = mu2 * (4*beta2 - 3*beta1) / (10*beta2 - 12*beta1 - 18)
    b1 = sqrt(mu2) * sqrt(beta1) * (beta2 + 3) / (10*beta2 - 12*beta1 - 18)
    b2 = (2*beta2 - 3*beta1 - 6) / (10*beta2 - 12*beta1 - 18)

    if isapprox(beta1, 0.0, rtol=0.0, atol=sqrt(eps()))
      if isapprox(beta2, 3.0, rtol=0.0, atol=sqrt(eps())) #Type 0
        return Normal(mu1, sqrt(mu2))
      elseif beta2 < 3 #Type II
        a1 = sqrt(mu2) / 2 * (-sqrt(-16 * beta2 * (2 * beta2 - 6)) / (2 * beta2 - 6))    
        m1 = - 1.0 / (2 * b2)
        return PearsonII(1+m1, mu1 - a1, 2*a1)
      else #beta2 > 3 #Type VII
        r = 6 * (beta2 - 1) / (2 * beta2 - 6)
        a = sqrt(mu2 * (r-1))
        return PearsonVII(1 + r, mu1, a / sqrt(1 + r))
      end
    elseif !isapprox((2*beta2 - 3*beta1 - 6), 0.0, rtol=0.0, atol=sqrt(eps()))
      k = 0.25 * beta1 * (beta2 + 3)^2 / ((4*beta2 - 3*beta1) * (2*beta2 - 3*beta1 - 6))
      if k < 0 #Type I
        a1 = sqrt(mu2) / 2 * ((-sqrt(beta1) * (beta2 + 3) - sqrt(beta1 * (beta2 + 3) ^ 2 - 4 * (4*beta2 - 3*beta1) * (2*beta2-3*beta1-6))) / (2*beta2-3*beta1-6))
        a2 = sqrt(mu2) / 2 * ((-sqrt(beta1) * (beta2 + 3) + sqrt(beta1 * (beta2 + 3) ^ 2 - 4 * (4*beta2 - 3*beta1) * (2*beta2-3*beta1-6))) / (2*beta2-3*beta1-6))
        if a1 > a2
          tmp = a1
          a1 = a2
          a2 = tmp
        end
        m1 =  -(sqrt(beta1) * (beta2+3) + a1*(10*beta2 - 12*beta1 - 18) / sqrt(mu2)) / 
        (sqrt(beta1 * (beta2 + 3)^2 -4*(4*beta2 - 3*beta1)*(2*beta2 - 3*beta1-6)))
        m2 = -(-sqrt(beta1) * (beta2+3) - a2*(10*beta2 - 12*beta1 - 18) / sqrt(mu2)) /
        (sqrt(beta1 * (beta2 + 3)^2 -4*(4*beta2 - 3*beta1)*(2*beta2 - 3*beta1-6)))
        lambda = mu1 - (a2-a1) * (1+m1) / (m1+m2+2)
        return PearsonI(1+m1, 1+m2, lambda, a2-a1)
      elseif isapprox(k, 1.0, rtol=0.0, atol=sqrt(eps())) #Type V
        return PearsonV(1/b2-1, mu1 - b1/(2*b2), -(b1-b1/(2*b2)) / b2) 
      elseif k > 1 #Type VI
        a1 = sqrt(mu2)/2*((-sqrt(beta1)*(beta2+3)-sqrt(beta1*(beta2+3)^2-4*(4*beta2-3*beta1)*
                                                       (2*beta2-3*beta1-6)))/(2*beta2-3*beta1-6))
        a2 = sqrt(mu2)/2*((-sqrt(beta1)*(beta2+3)+sqrt(beta1*(beta2+3)^2-4*(4*beta2-3*beta1)*
                                                       (2*beta2-3*beta1-6)))/(2*beta2-3*beta1-6))
        if a1 > a2
          tmp = a1
          a1 = a2
          a2 = tmp
        end
        m1 = (b1 + a1)/(b2*(a2-a1))
        m2 = - (b1 + a2)/(b2*(a2-a1))
        alpha = (a2-a1)
        lambda = mu1 + alpha * (m2+1)/(m2+m1+2)
        return PearsonVI(1+m2, -m2-m1-1, lambda, alpha)
      else # (k > 0) && (k < 1) #Type IV
        r = 6*(beta2-beta1-1)/(2*beta2-3*beta1-6)
        nu = -r*(r-2)*sqrt(beta1)/sqrt(16*(r-1)-beta1*(r-2)^2)
        alpha = sqrt(mu2*(16*(r-1)-beta1*(r-2)^2))/4
        lambda = mu1 - ((r-2)*sqrt(beta1)*sqrt(mu2))/4
        m = 1+r/2
        return PearsonIV(m, nu, lambda, alpha)
      end
    else # Type III
      m = b0 / (b1^2) - 1
      return PearsonIII(1+m, mu1 - b0/b1, b1)
    end
  end
end

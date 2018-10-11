using Plots

function pearson_plot(max_skewness::T=sqrt(14.0), max_kurtosis::T=24.0, n::Integer=301) where T <: Real
    t3fn(x) = 3.0 + 3.0/2.0 * x
    t5fn(x) = -3.0 * (16.0 + 13.0 * x + 2.0 * sqrt(64.0 + 48.0 * x + 12.0 * x^2.0 + x^3.0 )) / (x - 32.0)
    tnofn(x) = 1.0 + x;
    
    p = plot(xlim=(0, max_skewness^2), ylim=(0, max_kurtosis), xlab="beta_1", ylab="beta_2", leg=:bottomright)
    b1seq = range(0,stop=max_skewness^2, length=n)

    plot!(p, b1seq, tnofn.(b1seq), c=:black,  label="2-point-distr.", width=3)
    plot!(p, b1seq, t3fn.(b1seq), fillrange=[tnofn.(b1seq) t3fn.(b1seq)], c=:yellow, label="Pearson I")
    plot!(p, b1seq[b1seq.<32], fill(max_kurtosis, sum(b1seq.<32)), 
        fillrange=[t5fn.(b1seq[b1seq.<32]) fill(max_kurtosis, sum(b1seq.<32))], 
        c=:lightcyan, label="Pearson IV", width=2)
    plot!(p, b1seq[b1seq.<32], t5fn.(b1seq[b1seq.<32]),  
        fillrange=[t3fn.(b1seq[b1seq.<32]) t5fn.(b1seq[b1seq.<32])], c=:pink, label="Pearson VI")
    plot!(p, b1seq[b1seq.>=32], fill(max_kurtosis, sum(b1seq.>=32)),  
        fillrange=[t3fn.(b1seq[b1seq.>=32]) fill(max_kurtosis, sum(b1seq.>=32))], c=:pink, label="Pearson VI")

    plot!(p, b1seq[b1seq.<32], t5fn.(b1seq[b1seq.<32]), c=:blue, label="Pearson V", width=2)
    plot!(p, [fill(32, sum(b1seq .>= 32))], [fill(max_kurtosis, sum(b1seq.>=32))] , c=:blue, label="Pearson V", width=2)

    plot!(p, [0; 0], [1; 3], c=:orange, label="Pearson II", width=5)
    plot!(p, [0; 0], [3; max_kurtosis], c=:cyan, label="Pearson VII", width=5)
    plot!(p, b1seq, t3fn.(b1seq), c=:red, label="Pearson III", width=2)
    scatter!(p, [0],[3], c=:brown, label="Pearson 0")
    return p
end

pearson_plot(max_skewness::Real, max_kurtosis::Real, n::Integer) = 
    pearson_plot(promote(max_skewness, max_kurtosis)..., n)
pearson_plot(max_skewness::Integer, max_kurtosis::Integer, n::Integer) = 
    pearson_plot(Float64(max_skewness), Float64(max_kurtosis), n)

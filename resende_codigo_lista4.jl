using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve
using BenchmarkTools


@nbinclude("tauchen.ipynb") # para discretizar os choques
@nbinclude("ha_baseline.ipynb") # Elementos finitos + Galerkin
@nbinclude("ha_equilibrio.ipynb")


##########################################
## Distribuição invariante e Equilibrio ##
##########################################

@time R_C,inva_C, K_C = equilibrio_r(0.96,0.01,0.9,1.0001)[vec([1 2 6])]
@time R_D,inva_D, K_D = equilibrio_r(0.96,0.01,0.97,1.0001)[vec([1 2 6])]
@time R_E,inva_E, K_E = equilibrio_r(0.96,0.01,0.9,5)[vec([1 2 6])]
@time R_F,inva_F, K_F = equilibrio_r(0.96,0.05,0.9,1.0001)[vec([1 2 6])]


@time R_G,inva_G, K_G = equilibrio_r(0.96,0.05,0.9,5)[vec([1 2 6])]
@time R_H,inva_H, K_H = equilibrio_r(0.96,0.05,0.97,1.0001)[vec([1 2 6])]
@time R_I,inva_I, K_I = equilibrio_r(0.96,0.05,0.97,5)[vec([1 2 6])]
@time R_J,inva_J, K_J = equilibrio_r(0.96,0.01,0.97,5)[vec([1 2 6])]

S_C = exp.(disc_t(9; mu= 0, sigma = 0.01, rho= 0.9, m = 3)[1])
S_D = exp.(disc_t(9; mu= 0, sigma = 0.01, rho= 0.97, m = 3)[1])
S_E = exp.(disc_t(9; mu= 0, sigma = 0.01, rho= 0.9, m = 3)[1])
S_F = exp.(disc_t(9; mu= 0, sigma = 0.05, rho= 0.9, m = 3)[1])

100*R_F-100*R_C
A_C = LinRange(-minimum(S_C)/(1/0.96),4*maximum(S_C)/(1/0.96),1000)
A_D = LinRange(-minimum(S_D)/(1/0.96),4*maximum(S_D)/(1/0.96),1000)
A_E = LinRange(-minimum(S_E)/(1/0.96),4*maximum(S_E)/(1/0.96),1000)
A_F = LinRange(-minimum(S_F)/(1/0.96),4*maximum(S_F)/(1/0.96),1000)

plot(A_C,sum(inva_C,dims=1)', label = "ρ = 0.9, σ = 0.01, γ = 1.0001", title = "Aumento de persistência")
plot!(A_D,sum(inva_D,dims=1)', label = "ρ = 0.97, σ = 0.01, γ = 1.0001",ylims = (0,0.04))

plot(A_C,sum(inva_C,dims=1)', label = "ρ = 0.9, σ = 0.01, γ = 1.0001", title = "Aumento da Aversão ao risco")
plot!(A_E,sum(inva_E,dims=1)', label = "ρ = 0.9, σ = 0.01, γ = 5",ylims = (0,0.015))

plot(A_C,sum(inva_C,dims=1)', label = "ρ = 0.9, σ = 0.01, γ = 1.0001", title = "Aumento da variância")
plot!(A_F,sum(inva_F,dims=1)', label = "ρ = 0.9, σ = 0.05, γ = 1.0001",ylims = (0,0.04))


sum(inva_C,dims=1)[1]
sum(inva_D,dims=1)[1]
sum(inva_E,dims=1)[1]
sum(inva_F,dims=1)[1]

A_C
A_D
plot(S_C,sum(inva_C,dims=2), label = "ρ = 0.9, σ = 0.01, γ = 1.0001", title = "Distribuição marginal de renda")
plot!(S_D,sum(inva_D,dims=2), label = "ρ = 0.97, σ = 0.01, γ = 1.0001")
plot!(S_E,sum(inva_E,dims=2), label = "ρ = 0.9, σ = 0.01, γ = 5")
plot!(S_F,sum(inva_F,dims=2), label = "ρ = 0.9, σ = 0.05, γ = 1.0001")


plot(S_C,ones(9)*0.25, label = "ρ = 0.9, σ = 0.01, γ = 1.0001", title = "Espaço de choques para os parâmetros")
plot!(S_D,ones(9)*0.5, label = "ρ = 0.97, σ = 0.01, γ = 1.0001")
plot!(S_E,ones(9)*0.75, label = "ρ = 0.9, σ = 0.01, γ = 5")
plot!(S_F,ones(9), label = "ρ = 0.9, σ = 0.05, γ = 1.0001")
vline!([S_C[1] S_C[9]], label = nothing)
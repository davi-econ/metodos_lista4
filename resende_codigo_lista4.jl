using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve
using BenchmarkTools


@nbinclude("tauchen.ipynb") # para discretizar os choques
@nbinclude("poli_cheb.ipynb") # Funções para o polinomio de Chebyshev
@nbinclude("phi_consumo.ipynb") # Função de phi e consumo
@nbinclude("colocacao.ipynb") # Elementos finitos + Colocação
@nbinclude("ha_baseline.ipynb") # Elementos finitos + Galerkin

# parametros do modelo
β = beta = 0.96
γ = gamma = 1.0001
# r will clear the market
ρ = rho = 0.9
α = alpha = 1/3
σ = sigma = 0.01
δ = delta = 0.012


#####################
## grid de choques ##
#####################
Z, P = disc_t(9;mu= 0, sigma = 0.01, rho= 0.9, m = 3)
S = exp.(Z)

#############################
## discretização de ativos ##
#############################
grid_a = LinRange(-1,4,200)


########################
## problema do agente ##
########################
chute = zeros(length(grid_a),length(S))
@time v_ha, k_ha, c_ha = ha_baseline(grid_a;grid_z = S,P = P, chute = chute)

surface(v_ha)
surface(k_ha)
plot(c_ha)
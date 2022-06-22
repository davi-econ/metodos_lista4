using LinearAlgebra
using Distributions
using NBInclude
using Plots
using NLsolve
using BenchmarkTools


@nbinclude("tauchen.ipynb") # para discretizar os choques
@nbinclude("phi_consumo.ipynb") # Função de phi e consumo
@nbinclude("colocacao.ipynb") # Elementos finitos + Colocação
@nbinclude("ha_baseline.ipynb") # Elementos finitos + Galerkin
@nbinclude("ha_equilibrio.ipynb")
# parametros do modelo
β = beta = 0.96
γ = gamma = 3
# r will clear the market
ρ = rho = 0.9
α = alpha = 1/3
σ = sigma = 0.01
δ = delta = (1-β)/β
r = (1/β - 1)
r_0 = 0.04
sqrt(0.01*(1-0.3^2))
#####################
## grid de choques ##
#####################
Z, P = disc_t(9;mu= 0, sigma = 0.01, rho= 0.97, m = 3)
S = exp.(Z)

#############################
## discretização de ativos ##
#############################
ϕ = -minimum(S)/((1/β) - 1)
a_max = maximum(S)/(1/β - 1)
grid_a = LinRange(ϕ,a_max,1000)
########################
## problema do agente ##
########################

#############################
## Distribuição invariante ##
#############################
@time Q, λ =  dist_invariante(P; pk = k_ha)



surface(c_ha)
K,dist_a,dist_z, dist_m = hugget_capital_demand(λ;grid_a)
surface(dist_m)
plot(dist_a)
plot(dist_z)
erro2 =2
r_max = 0.0416666
r_min = 0

@time v_ha, k_ha, c_ha = ha_mg_acelerado([50 100 200 1000];grid_z = S,P,β = 0.96,γ = 1.0001,r = r_0)

plot(c_ha)

R_C = equilibrio_r(0.96,0.01,0.9,1.0001)[1]
R_D = equilibrio_r(0.96,0.01,0.97,1.0001)[1]
R_E = equilibrio_r(0.96,0.01,0.9,5)[1]
R_F = equilibrio_r(0.96,0.05,0.9,1.0001)[1]

R_D - R_C
K_new, dist_a, dist_z, matriz_M = hugget_capital_demand(λ;grid_a)
r_0
sum(matriz_M)
plot(dist_z)
surface(matriz_M)
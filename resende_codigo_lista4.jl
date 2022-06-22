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
Z, P = disc_t(9;mu= 0, sigma = 0.01, rho= 0.9, m = 3)
disc_t(9;mu= 0, sigma = 6, rho= 0.9, m = 3)[2]
S = exp.(Z)

#############################
## discretização de ativos ##
#############################
ϕ = -minimum(S)/(1+(1/β - 1))
a_max = minimum(S)/(1/β - 1)
grid_a = LinRange(ϕ,a_max,200)

########################
## problema do agente ##
########################
chute = zeros(length(grid_a),length(S))
@time v_ha, k_ha, c_ha = ha_baseline(grid_a;grid_z = S,P = P, chute = chute,α = 1/3,γ = 1.0001,β = 0.96, max_iter = 1000, r = 0.04)

surface(v_ha)
surface(k_ha)
surface(c_ha)
plot(c_ha)
Z
S
P
#############################
## Distribuição invariante ##
#############################
v_ha, k_ha, c_ha = ha_baseline(grid_a;grid_z = S,P = P, chute = chute,α = 1/3,γ = 3,β = 0.96, max_iter = 1000, r = 0.0384)

function dist_invariante(P;pk)
    nx = length(P[1,:])*length(pk[:,1])
    nz = length(P[1,:])
    Q = zeros(nx,nx)
    uns  = ones(nx,nx)
    um = ones(nx,1)
    Id1 = Matrix(I,nx,nx)
    for i in 1:nx
        for j in 1:nx
            a_i = div(j-1,nz) + 1 # índice de capital, para averiguar se é igrual à função política para aquele par
            b_i = (j - 1) % nz  + 1
            hoje_s = (i - 1)  % nz  + 1
            hoje_a= div(i-1,nz) + 1
            if a_i == pk[hoje_a,hoje_s]
                Q[i,j] = P[hoje_s,b_i]
            end
        end
    end
    λ = (Id1 - Q + uns)'\um
    return Q, λ
end

@time Q, λ =  dist_invariante(P; pk = k_ha)
plot(λ)
function hugget_capital_demand(dist_inva;grid_a)
    na = length(grid_a)
    nz = Int(length(dist_inva)/na)
    opa = reshape(dist_inva, nz,na)
    dist_a = sum(opa,dims = 1)
    dist_z = sum(opa,dims = 2)
    K = dot(grid_a,dist_a)
    return K, dist_a', dist_z
end

surface(c_ha)
K,dist_a,dist_z = hugget_capital_demand(λ;grid_a)
plot(dist_a)
plot(dist_z)
erro2 =2
r_max = 0.0416666
r_min = 0
r_0 = (r_max +r_min)/2
chute = ones(200,9)
sum(λ)
Q, λ =  dist_invariante(P; pk = k_ha)
K_new, dist_a, dist_z = hugget_capital_demand(λ;grid_a)
plot(dist_a)
grid_a = LinRange(-1,4,200)
K_0 = 0
while erro2 > 1e-7    
    v_ha, k_ha, c_ha = ha_baseline(grid_a;grid_z = S,P = P, chute = chute,α = 1/3,γ = 1.0001,β = 0.96, max_iter = 1000, r = r_0)
    Q, λ =  dist_invariante(P; pk = k_ha)
    K_new, dist_a, dist_z = hugget_capital_demand(λ;grid_a)
    r_new = α*real((K_new+ 0im)^(1-α)) + δ
    if K_new > 1e-5
        r_max = r_0
        r_min = r_min
        r_0 = (r_max +r_min)/2
        println("capital $K_new juros $r_0")
        erro2 = abs(r_new - r_0)
    elseif K_new < -1e-5
        r_max = r_max
        r_min = r_0
        r_0 = (r_max + r_min)/2
        println("capital $K_new juros $r_0")
        erro2 = abs(r_new - r_0)
    else
        println("capital $K_new juros $r_0")
        erro2 = abs(r_new - r_0)
    end
end


plot(dist_z)
K_new
using BlockMatrixHierarchy
using MosekTools
using JuMP
using MathOptInterface
using LinearAlgebra
using SparseArrays
using Ket
using Hypatia


d = 3
dAa = d+1
loc_d = d
num_parties::Int8  = 1
num_states::Int8  = 9
num_states_parties::Int8 = -1
inputs::Int8  = 1
outputs::Int8  = 9
level::Int8 = 1
state_kind::Int8 = 1
element_kind::Int8 = 1

arr=[]

eq = ["ρ1-ρ1"]
ineq = ["ρ1-ρ1"]

obj="ρ1-ρ1"

"""
z=exp(im*pi/4)
φ=(1+sqrt(5))/2
sic_fid = [sqrt(2), z*(1-z)*(φ^(3/2)+z'), (2-sqrt(2))*im, (z*(1-z)*(φ^(3/2)-z'))]
sic_norm = tr(ketbra(sic_fid))
targets = []
for a0 in [0,1,2,3]
    for a1 in [0,1,2,3]
        local tmp_state = shift(4)^a0*clock(4)^a1*sic_fid
        push!(targets, tmp_state*tmp_state'/sic_norm)
    end
end
"""

ω=exp(2*pi*im/d)
targets = Hermitian{ComplexF64, Matrix{ComplexF64}}[]
push!(targets, ketbra([0, 1, -1])/2)
push!(targets, ketbra([0, ω, ω^2])/2)
push!(targets, ketbra([0, ω^2, ω])/2)
push!(targets, ketbra([-1, 0, 1])/2)
push!(targets, ketbra([-ω^2, 0, ω])/2)
push!(targets, ketbra([-ω, 0, ω^2])/2)
push!(targets, ketbra([1, -1, 0])/2)
push!(targets, ketbra([ω, -ω^2, 0])/2)
push!(targets, ketbra([ω^2, -ω, 0])/2)

for i = 1:num_states
    targets[i]=embed(targets[i], dAa)
end

function upper_triangular(A::AbstractMatrix{String}; fill = "")
    mask = triu(trues(size(A)...))      # Bool mask: true on upper triangle
    return ifelse.(mask, A, fill)       # keep A where mask=true, else `fill`
end

custom_monomials="A11ρ1+A21ρ2+A31ρ3+A41ρ4+A51ρ5+A61ρ6+A71ρ7+A81ρ8+A91ρ9+ρ1A11+ρ2A21+ρ3A31+ρ4A41+ρ5A51+ρ6A61+ρ7A71+ρ8A81+ρ9A91"
custom_monomials=BlockMatrixHierarchy.parse_objective_terms(custom_monomials, 1, 1)
arr=[]
for n_loop=9:9
    ϵ = 0.01*(n_loop-1)

    global model, var_dict, variables, G, gamma = BlockMatSDP(dAa, dAa, num_parties, num_states, num_states_parties, Int8(0), Int8(0), level, eq, ineq, obj, state_kind, element_kind; Θ=:DirectSum, custom_monomials=custom_monomials)
    global SymG = symbolic_matrix(gamma)
    global SymG_up = upper_triangular(SymG; fill = "")
    global BigI = G[1:dAa, 1:dAa]


    for k=1:num_states
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ$(k)", SymG)]])) == 1)
    end

    for i=1:outputs
        for j=1:inputs
            for k=1:num_states
                @constraint(model, Hermitian(-var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)$(j)ρ$(i)", SymG)]] + var_dict[gamma[findfirst(x->x=="ρ$(i)", SymG)]]) in HermitianPSDCone() )
            end
        end
    end

    for k=1:num_states
        for j=1:inputs
            for i=1:outputs
                @constraint(model, 0 <= real(tr(var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)$(j)", SymG)]])) <= 1)
                @constraint(model, var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)$(j)", SymG)]] == var_dict[gamma[findfirst(x->x=="A$(i)$(j)ρ$(i)", SymG)]]')
            end
        end
    end

    for i=1:num_states
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ$(i)", SymG)]]*targets[i])) >= 1-ϵ )
    end

    model = trace_constraints(model, variables, var_dict)

    @objective(model, Max, real(tr(sum(var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)1", SymG)]] for i=1:outputs)))/num_states)
    optimize!(model)
    push!(arr, [objective_value(model), ϵ])
end




arr=[[el[1], el[2]] for el in arr]


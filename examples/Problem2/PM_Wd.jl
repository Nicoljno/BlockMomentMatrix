using BlockMatrixHierarchy
using MosekTools
using JuMP
using LinearAlgebra
using Ket
using Hypatia

function computational_basis_ρ(d::Integer)
    ρ = Matrix{ComplexF64}(I, d, d)               # identity
    return [ ρ[:,k] * ρ[k,:]'  for k in 1:d ]     # |k⟩⟨k|
end

function qft(d::Integer)
    ω = 2π * im / d
    F = [exp(ω*(j-1)*(k-1)) for j in 1:d, k in 1:d]
    return F / √d
end


d = 4
dAa = d+1
loc_d = d
num_parties::Int8  = 2
num_states::Int8  = 2*d #2*d
num_states_parties::Int8 = -1
inputs::Int8  = 1
outputs::Int8  = d #d
level::Int8 = 1
state_kind::Int8 = -1
element_kind::Int8 = 1

arr=[]

eq = ["ρ1-ρ1"]
ineq = ["ρ1-ρ1"]
#obj="1/3*(ρ1A11+ρ2A21+ρ3A31+ρ4B11+ρ5B21+ρ6B31)"

# d=2
#obj="(ρ1A11+ρ2A21+ρ3B11+ρ4B21)/2"
# d=3
#obj="(ρ1A11+ρ2A21+ρ3A31+ρ4B11+ρ5B21+ρ6B31)/3"
#obj="(ρ1A11+ρ2A21-ρ3A11-ρ3A21+ρ4B11+ρ5B21-ρ6B11-ρ6B21)"
# d=4
obj="(ρ1A11+ρ2A21+ρ3A31+ρ4A41+ρ5B11+ρ6B21+ρ7B31+ρ8B41)/4"
# d=5
#obj="(ρ1A11+ρ2A21+ρ3A31+ρ4A41+ρ5A51+ρ6B11+ρ7B21+ρ8B31+ρ9B41+ρ10B51)/5"
# d=6
#obj="(ρ1A11+ρ2A21+ρ3A31+ρ4A41+ρ5A51+ρ6A61+ρ7B11+ρ8B21+ρ9B31+ρ10B41+ρ11B51+ρ12B61)/6"

#ineq = ["ρ1A11", "ρ1B11", "ρ2A11", "ρ2B11", "ρ3A11", "ρ3B11", "ρ4A11", "ρ4B11"]
#ineq = [ineq; ["I/$d-ρ1A11", "I/$d-ρ1B11", "I/$d-ρ2A11", "I/$d-ρ2B11", "I/$d-ρ3A11", "I/$d-ρ3B11", "I/$d-ρ4A11", "I/$d-ρ4B11"] ]
#obj="1/8*(ρ1A11+ρ3A11+ρ1B11+ρ2B11-ρ4A11-ρ4B11-ρ2A11-ρ3B11)"

F = qft(d) #d
targets = computational_basis_ρ(d) #d
for i = 1:d #d
    push!(targets, F*targets[i]*F')
end
for i = 1:num_states
    targets[i]=embed(targets[i], dAa)
end

custom_monomials="A11ρ1+A21ρ2+A31ρ3+ρ1A11+ρ2A21+ρ3A31+B11ρ4+B21ρ5+B31ρ6+ρ4B11+ρ5B21+ρ6B31"
custom_monomials=obj
custom_monomials=BlockMatrixHierarchy.parse_objective_terms(custom_monomials, 1, -1)

function upper_triangular(A::AbstractMatrix{String}; fill = "")
    mask = triu(trues(size(A)...))      # Bool mask: true on upper triangle
    return ifelse.(mask, A, fill)       # keep A where mask=true, else `fill`
end

for n_loop=2:2
    ϵ = 0.01*(n_loop-1)

    global model, var_dict, variables, G, gamma = BlockMatSDP(dAa, dAa, num_parties, num_states, num_states_parties, inputs, outputs, level, eq, ineq, obj, state_kind, element_kind; Θ=:DirectSum, custom_monomials=custom_monomials)
    global SymG = symbolic_matrix(gamma)
    global SymG_up = upper_triangular(SymG; fill = "")
    global BigI = G[1:dAa, 1:dAa]

    for i=1:outputs
        for j=1:outputs
            if "A$(j)1B$(i)1" in SymG_up
                tmp = @variable(model, [1:dAa, 1:dAa] in HermitianPSDCone())
                @constraint(model, var_dict[gamma[findfirst(x->x=="A$(j)1B$(i)1", SymG)]] == tmp)
            end
        end
    end
    for j=1:outputs
        tmp = @variable(model, [1:dAa, 1:dAa] in HermitianPSDCone())
        @constraint(model, BigI - var_dict[gamma[findfirst(x->x=="A$(j)1", SymG)]] == tmp)
        tmp1 = @variable(model, [1:dAa, 1:dAa] in HermitianPSDCone())
        @constraint(model, BigI - var_dict[gamma[findfirst(x->x=="B$(j)1", SymG)]] == tmp1) 
    end
    for j=1:num_states
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]])) == 1)
        @constraint(model, sum(var_dict[gamma[findfirst(x->x=="ρ$(j)A$(i)1", SymG)]] for i=1:outputs)  == var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]])
        @constraint(model, sum(var_dict[gamma[findfirst(x->x=="ρ$(j)B$(i)1", SymG)]] for i=1:outputs)  == var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]])
    end
    for i=1:outputs
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)1ρ$(i)", SymG)]])  == tr(var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)1", SymG)]]))
        @constraint(model, tr(var_dict[gamma[findfirst(x->x=="ρ$(i+outputs)B$(i)1ρ$(i+outputs)", SymG)]])  == tr(var_dict[gamma[findfirst(x->x=="ρ$(i+outputs)B$(i)1", SymG)]]))
        @constraint(model, Hermitian(-var_dict[gamma[findfirst(x->x=="ρ$(i)A$(i)1ρ$(i)", SymG)]] + var_dict[gamma[findfirst(x->x=="ρ$(i)", SymG)]]) in HermitianPSDCone() )
        @constraint(model, Hermitian(-var_dict[gamma[findfirst(x->x=="ρ$(i+outputs)B$(i)1ρ$(i+outputs)", SymG)]] + var_dict[gamma[findfirst(x->x=="ρ$(i+outputs)", SymG)]]) in HermitianPSDCone() )
    end
    
    for j=1:num_states
        @constraint(model, Hermitian(var_dict[gamma[findfirst(x->x=="ρ$(j)", SymG)]] - var_dict[gamma[findfirst(x->x=="ρ$(j)ρ$(j)", SymG)]]) in HermitianPSDCone()) 
        for i=1:outputs
            @constraint(model, 0 <= real(tr(var_dict[gamma[findfirst(x->x=="ρ$(j)A$(i)1", SymG)]])) <= 1)
            @constraint(model, 0 <= real(tr(var_dict[gamma[findfirst(x->x=="ρ$(j)B$(i)1", SymG)]])) <= 1)
        end
    end

    @constraint(model, sum(var_dict[gamma[findfirst(x->x=="A$(i)1", SymG)]] for i=1:outputs) == BigI)
    @constraint(model, sum(var_dict[gamma[findfirst(x->x=="B$(i)1", SymG)]] for i=1:outputs) == BigI)
    @constraint(model, sum(var_dict[gamma[findfirst(x->x=="A$(i)1B$(j)1", SymG)]] for i=1:outputs for j=1:outputs) == BigI)

    ω = exp(2*π*im/3)
    for i=1:num_states
        @constraint(model, real(tr(var_dict[gamma[findfirst(x->x=="ρ$(i)", SymG)]]*targets[i])) >= 1-ϵ )
    end

    model = trace_constraints(model, variables, var_dict)

    #set_optimizer(model, Hypatia.Optimizer)
    #set_optimizer_attribute(model, "iter_limit", 150)
    optimize!(model)
    #push!(arr, (ε, value(obj)))
    push!(arr, [objective_value(model), ϵ])
end




arr=[[el[1], el[2]] for el in arr]
#FileIO.save("plot_quantum_RAC.jld2","arr",arr)
#FileIO.save("plot_classical_RAC_d$(d).jld2","arr",arr)


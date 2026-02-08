using BlockMatrixHierarchy
using MosekTools
using JuMP
using MathOptInterface
using LinearAlgebra
using SparseArrays
using Ket
using Hypatia

d=3
D=d^2
"""
z=exp(im*pi/4)
φ=(1+sqrt(5))/2
sic_fid = [sqrt(2), z*(1-z)*(φ^(3/2)+z'), (2-sqrt(2))*im, (z*(1-z)*(φ^(3/2)-z'))]
sic_norm = tr(ketbra(sic_fid))
targets = []
for a0 in [0,1,2,3]
    for a1 in [0,1,2,3]
        local tmp_state = shift(4)^a0*clock(4)^a1*sic_fid
        push!(targets, embed(tmp_state*tmp_state'/sic_norm, D))
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
for i = 1:d^2
    targets[i]=embed(targets[i], D)
end
P=random_povm(D,9)
arr=[]
for n_ϵ=1:100
    for n_loop=1:10
        global ϵ=0.0001*(n_ϵ-1)
        global model = Model(Mosek.Optimizer)

        targets_=[]
        for i=1:d^2
            push!(targets_, @variable(model, [1:D,1:D] in HermitianPSDCone()))
            @constraint(model, tr(targets_[i]) == 1)
            @constraint(model, real(tr(targets_[i]*targets[i])) >= 1-ϵ)
        end
        @objective(model, Max, real(tr(sum(P[i]*targets_[i] for i=1:d^2)/d^2)))
        optimize!(model)

        for i=1:d^2
            targets_[i]=value.(targets_[i])
        end

        global model = Model(Mosek.Optimizer)

        global P=[]
        for i=1:d^2
            push!(P, @variable(model, [1:D,1:D] in HermitianPSDCone()))
        end
        @constraint(model, sum(P[i] for i=1:d^2) == I)
        @objective(model, Max, real(tr(sum(P[i]*targets_[i] for i=1:d^2)/d^2)))
        optimize!(model)

        for i=1:d^2
            P[i]=value.(P[i])
        end
    end
    push!(arr, [objective_value(model), ϵ])
end
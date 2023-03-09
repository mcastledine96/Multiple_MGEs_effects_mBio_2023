cd(@__DIR__)
using Pkg
Pkg.activate(".")

using ClusterManagers, DifferentialEquations,CSV,DataFrames, LinearAlgebra, Distances,Distributions,Random, StatsBase
Random.seed!(parse(Int,ARGS[1]))
function LV(du,u,p,t)
    r, α,k = p
    n = length(u)
    @inbounds for i in 1:n
        du[i] = r[i]*u[i]
        @simd for j in 1:n
            du[i] -= r[i]*u[i]*α[j,i]*u[j]
        end
    end
end

function affect!(integrator)
    terminate!(integrator)
end
function condition(u,t,integrator)
  all(abs.(integrator.uprev .- u) .< 0.00001)
  end
cb =DiscreteCallback(condition,affect!)


nrep = 10000000
kendall =  Vector{Float64}(undef,nrep)
bray =  Vector{Float64}(undef,nrep)
alphas =  Vector{Float64}(undef,nrep)
stressors =  Vector{Float64}(undef,nrep)
ns =  Vector{Float64}(undef,nrep)
sv = Vector{Float64}(undef,nrep)
let model_batch = 0
    s = 0.1 + (parse(Int,ARGS[1]) ÷ 10)/10  
    for n in 2:2:20
        for alpha in 0.1:0.1:1     
            for i in 1:10000    
                    R =  abs.(randn(n))
                    A = rand(n,n) .*alpha/n
                    K = zeros(n)
                    A[diagind(A)] .= 1
                    P = (R,A,K)
                    u0 = abs.(randn(n))
                    prob = ODEProblem(LV,u0 ,(0,1000),P)
                    control = solve(prob,Tsit5(),callback=cb)[end]
                    control ./= sum(control)
                    Threads.@sync begin
                        for a in 1:1:10
                            Threads.@spawn begin
                                r = copy(R)
                                α = copy(A)
                                k = copy(K)
                                x = ones(n)
                                u = copy(u0)
                                for j in 1:a
                                    k .+= abs.(rand(Normal(0,s),n)) 
                                    push!(r,abs(randn()))
                                    push!(u,abs(randn()))
                                   
                                    α = vcat(α,zeros(1,size(α)[1]))
                                    α = hcat(α,zeros(size(α)[1]))
                                    α[end,:] .= abs.(rand(Normal(0,s),size(α)[2]))
                                    α[end,n+1:end] .= 1
                                    α[end,1:n] .= rand(truncated(Normal.(0,s),0,Inf),n)
                                    α[n+1:end,end] .= 1
                                    α[diagind(α)] .= 1
                                end
                                data_row = 10 * model_batch +a 
                                P = (r,α,k ./a)
                                prob = ODEProblem(LV,u,(0,1000),P)
                                sol = solve(prob,Tsit5(),callback=cb)
                                rep = sol[end][1:n] ./sum(sol[end][1:n])
                                kendall[data_row] = corkendall(control,rep)
                                bray[data_row] = evaluate(BrayCurtis(),control,rep)
                                ns[data_row] = n
                                sv[data_row] = s
                                stressors[data_row] = a
                                alphas[data_row] = alpha
                            end
                        end
                    end
                    model_batch+=1
            end
        end

df = DataFrame(hcat(kendall,
                    bray,
                    alphas,
                    stressors,
                    ns,
                    sv),[:kendall,:bray,:alpha,:stressors,:n,:stressvar])

CSV.write("SI_$(ARGS[1])_$(n).csv",df)
    end
end


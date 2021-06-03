module CoverFeatureMap

# set GLPK libs environment variables

using Random, GLPK, JuMP, LinearAlgebra

export cnp

function generate2plet(N, rho)
    x1 = vcat(zeros(N-1),1)
    x2 = vcat((randn(N-1) |> u-> (u/norm(u))*sqrt(1-rho^2)), rho)
    randmat = randn(N,N)
    Q, _ = qr(randmat)
    Q *= rand([1,-1]) * Matrix(I, (N, N))
    Q*x1, Q*x2, rand([-1,1])
end

function quadraticMap(x, linear::Bool, homogeneous::Bool, beta::Float64)
    N = length(x)
    phix = []
    for i in 1:N, j in i+1:N
        push!(phix, beta*x[i]*x[j])
    end
    if linear
        phix = vcat(phix, x)
    end
    if !homogeneous
        for i in 1:N
            push!(phix, beta*x[i]^2)
        end
    end
    phix/norm(phix)
end

function cubicMap(x, linear::Bool, homogeneous::Bool, beta::Float64)
    N = length(x)
    phix = []
    for i in 1:N, j in i+1:N, k in j+1:N
        push!(phix, beta*x[i]*x[j]*x[k])
    end
    if linear
        phix = vcat(phix, x)
    end
    if !homogeneous
        for i in 1:N
            push!(phix, beta*x[i]^3)
        end
    end
    phix/norm(phix)
end

# N : latent space dimension
# P : # of k-plets
# rho : overlap of k-plets
# nTrials : number of trials
function cnp(N::Int, P::Int, rho::Float64, 
            linear::Bool, homogeneous::Bool, beta::Float64, degree::Int,
            nTrials::Int;
            seed::Int = -1, precision::Float64=1e-6)

    seed > 0 && Random.seed!(seed)

    if degree == 2
        featureMap = x -> quadraticMap(x,linear,homogeneous,beta)
    elseif degree == 3
        featureMap = x -> cubicMap(x,linear,homogeneous,beta)
    else
        @error("Wrong degree")
    end
        
    #compute feature space dimension
    D = featureMap(ones(N)) |> length
    # println(D)


    overlaps = 0.
    sat = 0.
    for i in 1:nTrials
        # define linear program
        m = Model(GLPK.Optimizer)
    
        # set up model
        @variable(m, w[1:D])
        for j in 1:P
            x1, x2, sg = generate2plet(N, rho)
            y1 = featureMap(x1)
            y2 = featureMap(x2)
            overlaps += dot(y1,y2)

            @constraint(m, sg*dot(w,y1) >= precision)
            @constraint(m, sg*dot(w,y2) >= precision)
        end


        optimize!(m)
        # println(m)
        # println(JuMP.all_variables(m))
        if termination_status(m) == MOI.OPTIMAL
            sat += 1 
        end
    end
    
    overlap = overlaps / P / nTrials
    return N, P, rho, linear, homogeneous, beta, degree, nTrials, seed, precision, sat/nTrials, overlap
end

end # module

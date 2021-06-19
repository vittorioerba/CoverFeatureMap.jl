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

function generate3plet(N::Int64,  ρ12::Float64 , ρ23::Float64,  ρ13::Float64
                                                        ;seed::Int=-1)

    seed > 0 && Random.seed!(seed)

    ξ1 = vcat(zeros(N-1),1)
    ξ2 = vcat(zeros(N-2),[sqrt(1-ρ12^2), ρ12])

    temp1 =  (ρ23 -  ρ12* ρ13)/(sqrt(1 -  ρ12^2))
    temp2 = sqrt(1-  ρ13^2 - (( ρ23 -  ρ12*ρ13)^2/(1-  ρ12^2)))
    
    ξ3 = vcat(zeros(N-3),[temp2, temp1 , ρ13])
    
    randmat = randn(N,N)
    Q,R = qr(randmat)
    Q *= rand([1,-1])*Matrix(I, (N, N))
    
    return Q*ξ1, Q*ξ2, Q*ξ3, rand([-1,1])
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

function cnp_triplets(N::Int, P::Int, rho12::Float64, rho23::Float64, rho13::Float64, 
        linear::Bool, homogeneous::Bool, beta::Float64, degree::Int,
        nTrials::Int;
        seed::Int = -1, precision::Float64=1e-6)

    seed > 0 && Random.seed!(seed)

    if degree == 2
        featureMap = x -> quadraticMap(x, linear, homogeneous, beta)
    elseif degree == 3
        featureMap = x -> cubicMap(x, linear, homogeneous, beta)
    else
        @error("Wrong degree")
    end

    #compute feature space dimension
    D = featureMap(ones(N)) |> length
    # println(D)


    overlaps = [0., 0., 0.]
    sat = 0.
    for i in 1:nTrials
        # define linear program
        m = Model(GLPK.Optimizer)

        # set up model
        @variable(m, w[1:D])
        for j in 1:P
            x1, x2, x3, sg = generate3plet(N, rho12, rho23, rho13)
            y1 = featureMap(x1)
            y2 = featureMap(x2)
            y3 = featureMap(x2)
            overlaps[1] += dot(y1,y2)
            overlaps[2] += dot(y2,y3)
            overlaps[3] += dot(y1,y3)

            @constraint(m, sg*dot(w,y1) >= precision)
            @constraint(m, sg*dot(w,y2) >= precision)
            @constraint(m, sg*dot(w,y3) >= precision)
        end


        optimize!(m)
        # println(m)
        # println(JuMP.all_variables(m))
        if termination_status(m) == MOI.OPTIMAL
            sat += 1 
        end
    end

    overlap = [overlaps[i] / P / nTrials for i=1:3]
    return N, P, rho12, rho23, rho13, linear, homogeneous, beta, degree, nTrials, seed, precision, sat/nTrials, overlap[1], overlap[2], overlap[3]
end

end # module

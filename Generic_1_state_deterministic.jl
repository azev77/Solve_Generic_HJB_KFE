# generic solve_hjb( ) Implicit method. 
# TODO: ns states/nc controls 
# TODO: add uncertainty, nx shocks... Diffusions/Poisson/Jump-Diffusion

using LinearAlgebra, SparseArrays

# NGM ###############################
if 1==1 
    σ= 2.0; ρ = 0.05; δ = 0.05; A = 1.0; α= 0.3; 
    r(c; σ=σ)  = (c^(1-σ))/(1-σ)             # return fcn.      make it r(s,c)
    μ(s,c; α=α,δ=δ,A=A) = A*(s^α) - δ*s -c;  # Transition fcn.  ṡ=μ(s,c)
    dr(c;σ=σ) = c^(-σ)                       # deriv return fcn.
    FOC(v;σ=σ) = v^(-1.0/σ)                  # FOC, u'^-1()
    con0(s;α=α,δ=δ,A=A) = A*(s^α) - δ*s      # solve: ṡ=μ(s,c) for c(s) if ṡ=0. 

    s_ss = (α*A/(ρ+δ))^(1/(1-α))
    s_min = 0.001*s_ss
    s_max = 2.000*s_ss
    H = 10_000;
    s = LinRange(s_min, s_max, H)
    s = convert(Array, s)
    ds = (s_max-s_min)/(H-1)
end 

# con-sav ###############################
if 1==1 
    a_bar = -0.02;
    #a_bar = 0.1;
    σ= 2.0; ρ = 0.05; A=0.045; w = 0.1;
    r(c;σ=σ)  = (c^(1-σ))/(1-σ)            # return fcn.
    μ(s,c; α=α,δ=δ,A=A) = w + A*s - c;     # Transition fcn.  ṡ=μ(s,c)
    dr(c;σ=σ) = c^(-σ)                     # deriv return fcn.
    FOC(v;σ=σ) = v^(-1.0/σ)                # FOC, u'^-1()
    con0(s;α=α,δ=δ,A=A) = w + A*s          # solve: ṡ=μ(s,c) for c(s) if ṡ=0. 
    s_min = a_bar
    s_max = 1.0
    H = 10_000;
    s = LinRange(s_min, s_max, H)
    s = convert(Array, s)
    ds = (s_max-s_min)/(H-1)
end 

##################################
Δ = 1_000
maxit = 10_000
ε = 10e-6
dVf, dVb     = [zeros(H,1) for i in 1:2]
dV_Upwind, c = [zeros(H,1) for i in 1:2]

v0 = @. r(con0(s))/ρ #initial guess for V
v = v0
dist=[]
using Plots
for n=1:maxit
    println(n)
	V=v
    dV = (V[2:H]-V[1:H-1])/ds

    # forward difference
	dVf[1:H-1] = dV
	dVf[H]     = dr(con0(s_max))      # u'(c(s_max)) state constraint, for stability
    #dVf[H]= 0;

	# backward difference
	dVb[2:H] = dV
	dVb[1]   = dr(con0(s_min))        # state constraint, for stability
    #dVb[1] = (A*s_min)^(-σ)

    plot(dVf); plot!(dVb)

	I_concave = dVb .> dVf
    #scatter(I_concave) #1 everywhere EXCEPT @ last point H. 

    # consumption and savings with forward difference
    cf  = FOC.(dVf)
    μ_f = μ.(s, cf)

    # consumption and savings with backward difference
    cb  = FOC.(dVb)
    μ_b = μ.(s, cb)

    # c if k̇=0 & V = u'(c)
	c0  = con0.(s) 
    dV0 = dr.(c0)

    # Now to make a choice between forward and backward difference
    If = μ_f .> 0
    Ib = μ_b .< 0
    I0 = (1.0 .- If - Ib)

    dV_Upwind= dVf.*If + dVb.*Ib + dV0.*I0   
    c = FOC.(dV_Upwind)                       
    u = r.(c)

    # create the transition matrix
    X = -min.(μ_b,0)/ds
    Y = -max.(μ_f,0)/ds + min.(μ_b,0)/ds
    Z = max.(μ_f,0)/ds

    a1 = sparse(Diagonal((Y[:])))
    a2 = [zeros(1,H); sparse(Diagonal((X[2:H]))) zeros(H-1,1)]
    a3 = [zeros(H-1,1) sparse(Diagonal((Z[1:H-1]))); zeros(1,H)]
    AA = a1 + a2 + a3

    B = (ρ + 1/Δ)*sparse(I,H,H) - AA
    b = u + V./Δ

    # Solve: B V = b 
    # => (ρ + 1/Δ)*sparse(I,H,H)*V - AA*V = u + V./Δ
    # => (ρ + 1/Δ)*V = u + V./Δ +  AA*V
    # => (ρ)*V + V./Δ = u + V./Δ +  AA*V
    # => (ρ)*V = u  +  AA*V

    V = B \ b
    V_change = V-v
    v = V 

	push!(dist,findmax(abs.(V_change))[1])
	if dist[n] .< ε
		println("Value Function Converged Iteration=")
		println(n)
		break
	end
end

s_dot = μ.(s,c)
v_err = r.(c) + dV_Upwind.*s_dot - ρ.*v # approx @ borrowing constraint
# TODO: SIMULATE!

using Plots 
plot(dist, 
        #grid=false,
		xlabel="Iteration", ylabel="||V^{n+1} - V^n||",
		#ylims=(-0.001,0.030),
		legend=false, title="")
#png("Convergence")

plot(s, v_err, 
        #grid=false,
		xlabel="s", ylabel="Error in the HJB equation",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("HJB_error")

plot(s, v, 
        #grid=false,
		xlabel="k", ylabel="V(k)",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("Value_function_vs_k")

plot(s, c, 
        #grid=false,
		xlabel="k", ylabel="c(k)",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("c(k)_vs_k")

plot(s, s_dot, 
        #grid=false,
		xlabel="k", ylabel="s(k)",
		xlims=(s_min,s_max), title="", label="s(k)", legend=false)
plot!(s, zeros(H,1), label="", line=:dash)
#png("stateconstraint")

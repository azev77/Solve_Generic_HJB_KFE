# generic solve_hjb( ) Implicit method. 
# TODO: ns states/nc controls 
# TODO: add uncertainty, nx shocks...

using LinearAlgebra, SparseArrays, Plots

# NGM ###############################
if 1==1 
    σ= 2.0; ρ = 0.05; δ = 0.05; A = 1.0; α= 0.3; 

    # Closed form. Doesn't work.
    # σ= 4.0;  δ = 0.05; A = 1.0; α= 0.3; 
    # ρ = δ*(α*σ - 1) # restriction for closed form 
    # FOC(v;σ=σ) = abs(v)^(-1.0/σ)                  # FOC, u'^-1()

    r(c; σ=σ)  = (c^(1-σ))/(1-σ)             # return fcn.      make it r(s,c)
    μ(s,c; α=α,δ=δ,A=A) = A*(s^α) - δ*s -c;  # Transition fcn.  ṡ=μ(s,c)

    dr(c;σ=σ) = c^(-σ)                       # derivative of return fcn.
    FOC(v;σ=σ) = v^(-1.0/σ)                  # FOC, u'^-1()

    μ_inv(s,ṡ; α=α,δ=δ,A=A) = A*(s^α) - δ*s -ṡ; # Inv Trans c=μ(s,ṡ), solve ṡ=μ(s,c) for c
    # con0(s;α=α,δ=δ,A=A)     = A*(s^α) - δ*s      # solve: ṡ=μ(s,c) for c(s) if ṡ=0. 
    # μ_inv(s,0) = con0(s)

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
    μ(s,c; w=w,A=A) = w + A*s - c;         # Transition fcn.  ṡ=μ(s,c)

    dr(c;σ=σ) = c^(-σ)                     # deriv return fcn.
    FOC(v;σ=σ) = v^(-1.0/σ)                # FOC, u'^-1()

    μ_inv(s,ṡ; w=w,A=A) = w + A*s -ṡ;  # Inv Trans c=μ(s,ṡ), solve ṡ=μ(s,c) for c
    # μ_inv(s,0) = con0(s)
    #con0(s;α=α,δ=δ,A=A) = w + A*s          # solve: ṡ=μ(s,c) for c(s) if ṡ=0. 

    s_min = a_bar
    μ_inv(s_min,0) > 0 
    
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
v0 = @. r(μ_inv(s,0))/ρ #initial guess for V
v = v0
dist=[]

for n=1:maxit
    #println(n)
	V=v
    dV = (V[2:H]-V[1:H-1])/ds

    # forward difference
	dVf[1:H-1] = dV
	dVf[H]     = dr(μ_inv(s_max,0))      # u'(c(s_max)) state constraint, for stability
    #dVf[H]= 0;

	# backward difference
	dVb[2:H] = dV
	dVb[1]   = dr(μ_inv(s_min,0))        # state constraint, for stability
    #dVb[1] = (A*s_min)^(-σ)


	I_concave = dVb .> dVf
    # scatter(I_concave) #1 everywhere EXCEPT @ last point H. 

    # consumption and savings with forward difference
    cf  = FOC.(dVf)
    μ_f = μ.(s, cf)

    # consumption and savings with backward difference
    cb  = FOC.(dVb)
    μ_b = μ.(s, cb)

    # c if k̇=0 & V = u'(c)
	c0  = μ_inv.(s,0)   # con0.(s) 
    dV0 = dr.(c0)

    # Now to make a choice between forward and backward difference
    If = μ_f .> 0
    Ib = μ_b .< 0
    I0 = (1.0 .- If - Ib)

    dV_Upwind= dVf.*If + dVb.*Ib + dV0.*I0   
    c = FOC.(dV_Upwind)                       
    u = r.(c)

    # plot(dVf, lab="dVf") 
    # plot!(dVb, lab="dVb")
    # plot!(dV0, lab="dV0")
    # plot!(dV_Upwind, lab="dV_Upwind")

    # plot(cf, lab="dVf") 
    # plot!(cb, lab="dVb")
    # plot!(c0, lab="dV0")
    # plot!(c, lab="dV_Upwind")

    # plot(μ_f, lab="dVf") 
    # plot!(μ_b, lab="dVb")
    # plot!(zero(μ_f), lab="dV0")

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
    println(n, " ", dist[n])
	if dist[n] .< ε
		println("Value Function Converged Iteration=")
		println(n)
		break
	end
end

s_dot = μ.(s,c)
v_err = r.(c) + dV_Upwind.*s_dot - ρ.*v # approx @ borrowing constraint
# TODO: SIMULATE!

plot(dist, 
		xlabel="Iteration", ylabel="||V^{n+1} - V^n||",
		#ylims=(-0.001,0.030),
		legend=false, title="")
#png("Convergence")

plot(s, v_err, 
		xlabel="s", ylabel="Error in the HJB equation",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("HJB_error")

plot(s, v, 
		xlabel="k", ylabel="V(k)",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("Value_function_vs_k")

plot(s, c, 
		xlabel="k", ylabel="c(k)",
		xlims=(s_min,s_max),
		legend=false, title="")
#png("c(k)_vs_k")

plot(s, s_dot, 
		xlabel="k", ylabel="s(k)",
		xlims=(s_min,s_max), title="", label="s(k)", legend=false)
plot!(s, zeros(H,1), label="", line=:dash)
#png("stateconstraint")

# Simulate. 
# s_0 = 0.10;
# c_0 # = interpolate(s_0, (s,c(s)))
# sp = μ(s_0, c)






# Special case of NGM has closed-form sol. Compare!
σ= 4.0;  δ = 0.05; A = 1.0; α= 0.3; 
ρ = δ*(α*σ - 1) # restriction for closed form 

# Policy: c = (1-s)(k^α)
cp(s) = (1-((α*δ)/(ρ+δ)))*s^α
s_d(s) = (s^α) - δ*s - cp(s)
plot(s, cp)
plot!(s, s_d)

#
s_ss = (α*A/(ρ+δ))^(1/(1-α))
z_ss = (s_ss)^(1-α)
c_ss = (s_ss)^(α)  - δ*s_ss

s_0 = 5;
z_0 = (s_0)^(1-α)

λ=(1-α)*δ
z(t) = z_ss + exp(-λ*t)*(z_0 - z_ss)
kk(t) = z(t)^(1/(1-α))
cc(t) = (1-((α*δ)/(ρ+δ)))*kk(t)^α

t_sim = 0.01:0.01:100
plot(legend=:topleft)
plot!(t_sim, z, lab="z")
plot!(t_sim, kk, lab="k")
plot!(t_sim, cc, lab="c")
plot!([z_ss s_ss c_ss],  seriestype = :hline, lab="", color="grey")


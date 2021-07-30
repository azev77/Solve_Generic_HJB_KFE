# generic solve_hjb( ) Implicit method. 
# TODO: ns states/nc controls 
# TODO: add uncertainty, nx shocks...

# fcn NGM() returns struct NGM; 
# HJB_options() returns struct HJB_options;
# solve_HJB: (m, opt, v) -> dist, v_err, v, dV_Upwind, c_Upwind, μ_Upwind
# Plot error/value/policy/
# Simulate w/ interp.

# optimal stopping w/ LCPsolve.jl

using LinearAlgebra, SparseArrays, Plots

_m(x) = max(x,1e-10)
_m(x) = max(x,1e-8)

# NGM ###############################
struct NGM
    ρ; σ;            # preference param 
    δ; A; α;         # tech param 
    #
    r                # return fcn. r(s,c)
    μ                # transition fcn.  ṡ=μ(s,c)
    dr               # derivative of return fcn.
    FOC              # FOC: solve dr(c)=V_s. c=u'^-1(V_s)
    μ_inv            # inv transition c=μ_inv(s,ṡ), solve ṡ=μ(s,c) for c
    #
    s_min; s_max; H; s; ds;
    dVf; dVb; dV_Upwind; c_Upwind;
end

function NGM(; 
    ρ = 0.05, σ= 2.0,  
    δ = 0.05, A = 1.0, α= 0.3,
    r = c -> (c^(1-σ))/(1-σ),               # return fcn.      make it r(s,c)
    μ = (s,c) -> A*(s^α) -δ*s -c,           # Transition fcn.  ṡ=μ(s,c)
    dr = c -> c^(-σ),                       # derivative of return fcn.
    FOC = v -> (max(v,1e-10))^(-1.0/σ),     # FOC: solve dr(c)=V_s. c=u'^-1(V_s)
    μ_inv = (s,ṡ) -> A*(s^α) - δ*s -ṡ,      # Inv Trans c=μ_inv(s,ṡ), solve ṡ=μ(s,c) for c
    s_min = 0.001*(α*A/(ρ+δ))^(1/(1-α)),
    s_max = 2.000*(α*A/(ρ+δ))^(1/(1-α)),
    H = 10_000,
    s = collect(LinRange(s_min, s_max, H)),
    ds = (s_max-s_min)/(H-1), 
    dVf = zeros(H,1),
    dVb = zeros(H,1),
    dV_Upwind = zeros(H,1),
    c_Upwind = zeros(H,1),
    )
    return NGM(ρ,σ,δ,A,α,r,μ,dr,FOC,μ_inv,s_min,s_max,H,s,ds,dVf,dVb,dV_Upwind,c_Upwind)
end
mod = NGM()
mod = NGM(ρ = 0.05)
mod = NGM(A = 0.01)

# con-sav ###############################
struct CS
    ρ; σ;            # preference param 
    A; w;            # tech param 
    #
    r                # return fcn. r(s,c)
    μ                # transition fcn.  ṡ=μ(s,c)
    dr               # derivative of return fcn.
    FOC              # FOC: solve dr(c)=V_s. c=u'^-1(V_s)
    μ_inv            # inv transition c=μ_inv(s,ṡ), solve ṡ=μ(s,c) for c
    #
    s_min; s_max; H; s; ds;
    dVf; dVb; dV_Upwind; c_Upwind;
end

function CS(; 
    ρ = 0.05, σ= 2.0,  
    A = 0.045, w = 0.1,                     # A = interest rate 
    r = c -> (_m(c)^(1-σ))/(1-σ),           # return fcn.      make it r(s,c)
    μ = (s,c) -> w + A*s - c,               # Transition fcn.  ṡ=μ(s,c)
    dr = c -> _m(c)^(-σ),                       # derivative of return fcn.
    FOC = v -> (_m(v))^(-1.0/σ),     # FOC: solve dr(c)=V_s. c=u'^-1(V_s)
    μ_inv = (s,ṡ) -> w + A*s -ṡ,      # Inv Trans c=μ_inv(s,ṡ), solve ṡ=μ(s,c) for c
    s_min = -0.02,
    s_max = 1.0,
    H = 10_000,
    s = collect(LinRange(s_min, s_max, H)),
    ds = (s_max-s_min)/(H-1), 
    dVf = zeros(H,1),
    dVb = zeros(H,1),
    dV_Upwind = zeros(H,1),
    c_Upwind = zeros(H,1),
    )
    return CS(ρ,σ,A,w,r,μ,dr,FOC,μ_inv,s_min,s_max,H,s,ds,dVf,dVb,dV_Upwind,c_Upwind)
end

# NBC 
# μ_inv(a_bar,0) > 0 ⟺ w + A*a_bar >0 ⟺ a_bar > -w/A
#a_bar = -w/A # Natural Borrowing Constraint 
#a_bar = -w/A + eps()

mod = CS()
mod = CS(w = 0.1, s_min = -0.1/0.045 + eps())

# verify things are well defined @ corners 
mod.μ_inv(mod.s_min,0) > 0
mod.s_max, mod.μ_inv(mod.s_max,0), mod.dr(mod.μ_inv(mod.s_max,0))    
mod.s_min, mod.μ_inv(mod.s_min,0), mod.dr(mod.μ_inv(mod.s_min,0))    


##################################
# options 
struct HJB_options; Δ; maxit; ε; dist; end
HJB_options(;Δ=1_000,maxit=10_000,ε=1e-6,dist=[]) = HJB_options(Δ,maxit,ε,dist)
HJB_options()


##################################
# solve_HJB(m, opt, v) 
function solve_HJB(m, opt, v)
    # unpack 
    # σ,δ,A,α = m.σ, m.δ, m.A, m.α     # Do not need. 
    ρ = m.ρ
    r,μ,dr,FOC,μ_inv = m.r, m.μ, m.dr, m.FOC, m.μ_inv
    s_min,s_max,H,s,ds = m.s_min, m.s_max, m.H, m.s, m.ds
    dVf,dVb,dV_Upwind,c_Upwind = m.dVf, m.dVb, m.dV_Upwind, m.c_Upwind
    #
    Δ,maxit,ε,dist = opt.Δ, opt.maxit, opt.ε, opt.dist

    # VFI 
    for n=1:maxit
        #println(n)
        V=v
        dV = (V[2:end] - V[1:end-1])/ds  
    
        # dVf = [dV; dr(μ_inv(s_max,0))]
        # dVb = [dr(μ_inv(s_min,0)); dV]
    
        # forward difference: if ṡ>0
        dVf[1:end-1] = dV
        dVf[end]     = dr(μ_inv(s_max,0))    # u'(c(s_max)) state const, for stab #dVf[H]=0
        cf           = FOC.(dVf)             # choice with forward difference
        μ_f          = μ.(s, cf)             # ṡ      with forward difference
        If = μ_f .> 0
    
        # backward difference: if ṡ<0
        dVb[2:end] = dV
        dVb[1]     = dr(μ_inv(s_min,0))      # state const, for stab #dVb[1] = (A*s_min)^(-σ)
        cb  = FOC.(dVb)                      # choice with backward difference
        μ_b = μ.(s, cb)
        Ib = μ_b .< 0
    
        # neither difference: if ṡ=0
        c0  = μ_inv.(s,0)        # c if ṡ=0
        dV0 = dr.(c0)
        # μ_0 = μ.(s, c0)          # μ_0 == zero(s)
        I0  = (1.0 .- If - Ib)   # choice betw forward & backward difference
    
        # I_concave = dVb .> dVf
        # scatter(I_concave) #1 everywhere EXCEPT @ last point H. 
     
        dV_Upwind = dVf.*If + dVb.*Ib + dV0.*I0   
        c_Upwind  = FOC.(dV_Upwind)    
        μ_Upwind  = μ.(s, c_Upwind)                                            
        u = r.(c_Upwind)
    
        # plot(dVf, lab="dVf") 
        # plot!(dVb, lab="dVb")
        # plot!(dV0, lab="dV0")
        # plot!(dV_Upwind, lab="dV_Upwind")
    
        # plot(cf, lab="dVf") 
        # plot!(cb, lab="dVb")
        # plot!(c0, lab="dV0")
        # plot!(c_Upwind, lab="dV_Upwind")
    
        # plot(μ_f, lab="dVf") 
        # plot!(μ_b, lab="dVb")
        # plot!(zero(μ_f), lab="dV0")
        # plot!(μ_Upwind, lab="dV_Upwind")
    
    
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
    μ_Upwind  = μ.(s, c_Upwind)
    v_err = r.(c_Upwind) + dV_Upwind.*μ_Upwind - ρ.*v # approx @ borrowing constraint

    return dist, v_err, v, dV_Upwind, c_Upwind, μ_Upwind
end

##################################
# DO IT. 
mod = NGM()

# NGM Closed form. Need: ρ = δ*(α*σ - 1)
mod = NGM(σ= 4.0, δ = 0.05, A = 1.0, α= 0.3, ρ = 0.05*(0.3*4 - 1))

mod = CS()
# verify TVC: ρ-(1-σ)*A >0
# mod.ρ - (1-mod.σ)*mod.A   # AZ: always true if σ>1

mod.w, mod.s_min, -mod.w/mod.A
# AZ: lower s_min, the longer it takes to hit constraint. Bigger T_sim
mod = CS(s_min=-.5)
mod = CS(s_min=-1.5)
mod = CS(s_min=-2)
mod = CS(s_min=-0.1/0.045 + eps())
mod = CS(s_min=-0.1/0.045 -.25)

mod = CS(w=0.0, s_min=0.015) # binding
mod = CS(w=0.0, s_min=0.0)   # NBC
mod = CS(w=0.0, s_min=-0.01)  # loose 
#mod = CS(w=0.0, s_min=-0.10)  # loose. Inaccurate when s close to 0 

# AZ: got errors before bc w=0 & s_min =-.02

#
opt = HJB_options(ε=1e-7)
v0 = @. mod.r(mod.μ_inv(mod.s,0))/mod.ρ #v0 = @. r(μ_inv(s,0))/ρ #initial guess for V

dist, v_err, v, dV_Upwind, c_Upwind, μ_Upwind = solve_HJB(mod, opt, v0)

##################################
# Plot: dist, v_err, v, c_Upwind, μ_Upwind
plot(dist, 
		xlabel="Iteration", ylabel="||V^{n+1} - V^n||",
		legend=false, title="")
plot(mod.s, v_err, 
		xlabel="s", ylabel="Error in the HJB equation",
		legend=false, title="")
plot(mod.s, v, 
		xlabel="s", ylabel="V(s)",
		legend=false, title="")
plot(mod.s, c_Upwind, 
		xlabel="s", ylabel="c(s)",
		legend=false, title="")
plot(mod.s, μ_Upwind, 
		xlabel="s", ylabel="ṡ(s)",
		title="", label="ṡ(s)", legend=false);
plot!([0.0],  seriestype = :hline, lab="", color="grey", l=:dash)

##################################
# Simulate: s, c, ṡ
# problem: c(s) only on finite grid s => interpolate
using Interpolations
ĉ = LinearInterpolation(mod.s, c_Upwind[:], extrapolation_bc = Line())

Δt = 0.001; 
Δt = 0.01;
T_sim = 200; 
T_sim = 500;
T_sim = 1_000;

T_sim = 3_500;

t_sim = 0.0:Δt:T_sim
s_sim, c_sim, ṡ_sim = [zeros(length(t_sim),1) for i in 1:3]

s_0=mod.s[1] # start below ss 
s_0=0.5
#s_0=0.75

s_sim[1] = s_0 
for i in 2:length(t_sim)
    c_sim[i-1] = ĉ(s_sim[i-1]) 
    ṡ_sim[i-1] = mod.μ(s_sim[i-1], c_sim[i-1])
    s_sim[i]   = s_sim[i-1] + Δt * ṡ_sim[i-1]
    # s_sim[i] = s_sim[i-1] + Δt * μ(s_sim[i-1], c_sim[i-1])
end
#
ix = 1:(length(t_sim)-1)
plot();
plot!(t_sim[ix], s_sim[ix], lab = "s");
plot!(t_sim[ix], c_sim[ix], lab = "c");
plot!(t_sim[ix], ṡ_sim[ix], lab = "ṡ");
plot!([mod.s_min],  seriestype = :hline, lab="", color="grey", l=:dash)

# Verify: HACT Fig 2 Cara/F3 CRRA, hit borrowing constraint in finite time 



########################################################
# Special case of NGM has closed-form sol. Compare Policy & Sim 
# Policy: c = (1-s)(k^α)
cp(s) = (1-((α*δ)/(ρ+δ)))*s^α
s_d(s) = (s^α) - δ*s - cp(s)
plot(legend=:topleft)
plot!(s, cp,      lab="c: closed form")
plot!(s, c_Upwind,lab="c: Upwind")
plot!(s, s_d,     lab="ṡ: closed form")
plot!(s, s_dot,   lab="ṡ: Upwind")

#
s_ss = (α*A/(ρ+δ))^(1/(1-α))
z_ss = (s_ss)^(1-α)
c_ss = (s_ss)^(α)  - δ*s_ss

#s_0 = 5;
z_0 = (s_0)^(1-α)

λ=(1-α)*δ
z(t) = z_ss + exp(-λ*t)*(z_0 - z_ss)
kk(t) = z(t)^(1/(1-α))
cc(t) = (1-((α*δ)/(ρ+δ)))*kk(t)^α

#t_sim = 0.01:0.01:100
t_sim = time[ix]

plot(legend=:topleft)
#plot!(t_sim, z, lab="z")
plot!(t_sim, kk, lab="s: closed form")
plot!(t_sim, s_sim[ix], lab="s: Upwind")
plot!(t_sim, cc, lab="c: closed form")
plot!(t_sim, c_sim[ix], lab="c: Upwind")
#plot!([z_ss s_ss c_ss],  seriestype = :hline, lab="", color="grey")
plot!([s_ss c_ss],  seriestype = :hline, lab="", color="grey")



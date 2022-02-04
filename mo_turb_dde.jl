using DifferentialEquations
using Parameters
using CairoMakie # use Makie instead, well worth the effort to learn!!!!

# The basic turbidostat model with time-delay effects (for N_2)

@with_kw struct parameters
    c₁::Float64 = 0.13
    c₂::Float64 = 0.11
    c₃::Float64 = 0.05
    a::Float64 = 0.95
    b::Float64 = 7.5
    c::Float64 = 0.34
    d::Float64 = 8.2
    e::Float64 = 0.37
    f::Float64 = 1.52
    OD_chlamy::Float64 = 1.0
    OD_bacteria::Float64 = 1.0
    α::Float64 = 1.5
    g₁::Float64 = 10.0
    g₂::Float64 = 12.0
    d₁::Float64 = 0.005
    d₂::Float64 = 0.05
    pop_c::Float64 = 0.025
    pop_sc::Float64 = 0.025
    Dil_const::Float64 = 0.072
    N₁⁰::Float64 = 10.0
    N₂⁰::Float64 = 0.2
    C⁰::Float64 = 0.01
    SC⁰::Float64 = 0.1
    threshold::Float64 = 10.0
    Δ⁰::Int = 0
    tau::Int = 15
end

model_constants = parameters()   # Create the default parameters

# make functions one liners for readability
sub_asm(a, b, X) = a * X / (b + X)
sub_asm_lin(a, X) = a * X
sub_asm_cat(a, b, X) = a * (X^2) / (b^2 + X^2)
sub_asm_exp(a, b, X) = a * (1 - exp(-b * X))
growth_sat(a, X) = a / (a + X)
dilution(X, Y, a, b, d) = a * X + b * Y + d

out = zeros(5) # Define a cache variable
# Define the DDE system, parameter list p can be retrieved from the previous settings
function dde_sys(du, u, h, p, t)
    hist = h(p, t - p.tau)[3] # updates out to be the correct history function
    
    du[1] =
        u[5] * dilution(u[3], u[4], p.pop_c, p.pop_sc, 0) * (p.N₁⁰ - u[1]) -
        p.c₁ * sub_asm(p.a, p.b, u[1]) * u[3] / p.OD_chlamy -
        p.c₂ * sub_asm(p.c, p.d, u[1]) * u[3] / p.OD_bacteria

    du[2] =
        p.α * hist -
        p.c₃ * sub_asm(p.e, p.f, u[2]) * u[4] / p.OD_bacteria -
        u[5] * u[2] * dilution(u[3], u[4], p.pop_c, p.pop_sc, 0)

    du[3] =
        u[3] * (
            p.c₁ * sub_asm(p.a, p.b, u[1]) -
            u[5] * dilution(u[3], u[4], p.pop_c, p.pop_sc, 0) -
            p.d₁
        )

    du[4] =
        u[4] * (
            p.c₂ * sub_asm(p.c, p.d, u[1]) +
            p.c₃ * sub_asm(p.e, p.f, u[2] - p.d₂) -
            u[5] * dilution(u[3], u[4], p.pop_c, p.pop_sc, 0)
        )

    du[5] = 0
end

activate_dilution_condition(u, t, integrator) =
    integrator.u[5] < 1 && integrator.u[3] + integrator.u[4] > integrator.p.threshold

activate_dilution_affect!(integrator) = integrator.u[5] = 1

stop_dilution_condition(u, t, integrator) =
    integrator.u[5] > 0 && integrator.u[3] + integrator.u[4] < 0.9 * integrator.p.threshold

stop_dilution_affect!(integrator) = integrator.u[5] = 0.0

cb1 = DiscreteCallback(activate_dilution_condition, activate_dilution_affect!)
cb2 = DiscreteCallback(stop_dilution_condition, stop_dilution_affect!)
cbs = CallbackSet(cb1, cb2)

h(p, t) = ones(5)
tau = model_constants.tau
lags = [tau]

tspan = (0.0, 1000.0)
u0 = [
    model_constants.N₁⁰,
    model_constants.N₂⁰,
    model_constants.C⁰,
    model_constants.SC⁰,
    model_constants.Δ⁰,
]

prob = DDEProblem(dde_sys, u0, h, tspan, model_constants; constant_lags=lags)


alg = MethodOfSteps(Tsit5())
                                                                        
sol = solve(prob,alg, callback=cbs, saveat=0.1) 
                                                                        
plot_layout = @layout [
    grid(2,2)
    a{0.2h}
]

plot(sol,
    dpi = 500,
    size = (600,600),
    layout = plot_layout,
    title = ["N₁" "N₂" "Chlamy" "SynCom" "Dilution Switch"],
    label = ["" "" "" "" ""]
    )

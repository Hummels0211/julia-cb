using DifferentialEquations
using Parameters
using Plots

# The basic turbidostat model

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
end

model_constants = parameters()

# make functions one liners for readability
sub_asm(a, b, X) = a * X / (b + X)
sub_asm_lin(a, X) = a * X
sub_asm_cat(a, b, X) = a * (X^2) / (b^2 + X^2)
sub_asm_exp(a, b, X) = a * (1 - exp(-b * X))
growth_sat(a, X) = a / (a + X)
dilution(X, Y, a, b, d) = a * X + b * Y + d

# Define the ODE system, parameter list p can be retrieved from the previous settings
function ode_sys(du, u, p, t)
    du[1] =
        u[5] * dilution(u[3], u[4], p.pop_c, p.pop_sc, 0) * (p.N₁⁰ - u[1]) -
        p.c₁ * sub_asm(p.a, p.b, u[1]) * u[3] / p.OD_chlamy -
        p.c₂ * sub_asm(p.c, p.d, u[1]) * u[3] / p.OD_bacteria

    du[2] =
        p.α * u[3] -
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

# Define the stochasity adding to the ODE system
function σ_ode_sys(du, u, p, t)
    du[1] = 0.03*u[1]
    du[2] = 0.05*u[2]
    du[3] = 0.05*u[3]
    du[4] = 0.05*u[4]
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

tspan = (0.0, 1000.0)
u0 = [
    model_constants.N₁⁰,
    model_constants.N₂⁰,
    model_constants.C⁰,
    model_constants.SC⁰,
    model_constants.Δ⁰,
]

# Choose the random process for simulating the SDE system
W = WienerProcess(0.0,0.0,0.0)
# OU = OrnsteinUhlenbeckProcess(1.0, 1.2, 0.5, 0.0, 0.0)

prob = SDEProblem(ode_sys, σ_ode_sys, u0, tspan, model_constants, noise=W)
sol = solve(
    prob,
    SRIW1(),
    callback=cbs,
    saveat=0.1
)


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

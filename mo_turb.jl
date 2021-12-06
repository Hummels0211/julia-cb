using DifferentialEquations
using Plots

model_constants = ( # use a named tuple instead
    c₁ = 0.1,
    c₂ = 0.35,
    c₃ = 0.31,
    a = 2.5,
    b = 7.5,
    c = 0.34,
    d = 10.2,
    e = 0.37,
    f = 0.52,
    OD_chlamy = 1,
    OD_bacteria = 1,
    α = 1.5,
    g₁ = 10,
    g₂ = 12,
    d₁ = 0.005,
    d₂ = 0.05,
    pop_c = 0.06,
    pop_sc = 0.06,
    Dil_const = 0.072,
    N₁⁰ = 10.0,
    N₂⁰ = 0.2,
    C⁰ = 0.01,
    SC⁰ = 0.1,
    threshold = 10,
    Δ⁰ = 0,  # dilution state (0 inactive; 1, active)
)

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

activate_dilution_condition(u, t, integrator) =
    integrator.u[5] < 1 && integrator.u[3] + integrator.u[4] > integrator.p.threshold

activate_dilution_affect!(integrator) = integrator.u[5] = 1

stop_dilution_condition(u, t, integrator) =
    integrator.u[5] > 0 && integrator.u[3] + integrator.u[4] < 0.9 * integrator.p.threshold

stop_dilution_affect!(integrator) = integrator.u[5] = 0.0

cb1 = DiscreteCallback(activate_dilution_condition, activate_dilution_affect!)
cb2 = DiscreteCallback(stop_dilution_condition, stop_dilution_affect!)
cbs = CallbackSet(cb1, cb2)

tspan = (0.0, 100.0)
u0 = [
    model_constants.N₁⁰,
    model_constants.N₂⁰,
    model_constants.C⁰,
    model_constants.SC⁰,
    model_constants.Δ⁰,
]
prob = ODEProblem(ode_sys, u0, tspan, model_constants)
sol = solve(
    prob,
    callback=cbs,
    saveat=0.1
)


plot_layout = @layout [
    grid(2,2)
    a{0.2h}
]

plot(sol,
    dpi = 500,
    size = (600,400),
    layout = plot_layout,
    label = ["N₁" "N₂" "Chlamy" "SynCom" "Dilution Switch"],
    legend = :left)

# using DataFrames

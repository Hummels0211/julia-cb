using DifferentialEquations
using Parameters
using Plots

# The basic turbidostat model

@with_kw struct parameters
    c₁::Float64 = 0.13   # The assimilation rate of C
    c₂::Float64 = 0.11   # The assimilation rate of SC1
    c₃::Float64 = 0.05   # The non-metabolic assimilation rate of SC1
    c₄::FLoat64 = 0.12   # The assimilation rate of SC2
    c₅::Float64 = 0.06   # The non-metabolic assimilation rate of SC2
    a::Float64 = 0.95    # The assimilation constants of C
    b::Float64 = 7.5
    c_1::Float64 = 0.34   # The assimilation constants of SC1
    d_1::Float64 = 8.2
    e_1::Float64 = 0.37   # The non-metabolic assimilation constants of SC1
    f_1::Float64 = 1.52
    c_2::Float64 = 0.24   # The assimilation constants of SC2
    d_2::Float64 = 10.2
    e_2::Float64 = 0.37   # The non-metabolic assimilation constants of SC2  
    f_2::Float64 = 1.52
    OD_chlamy::Float64 = 1.0
    OD_bacteria::Float64 = 1.0
    α::Float64 = 1.5
    g₁::Float64 = 10.0
    g₂::Float64 = 12.0
    g₃::Float64 = 7.0
    d₁::Float64 = 0.005   # The death rate of Chlamy
    d₂::Float64 = 0.005   # The death rate of SC1
    d₃::Float64 = 0.005   # The death rate of SC2
    pop_c::Float64 = 0.025
    pop_sc_1::Float64 = 0.025
    pop_sc_2::Float64 = 0.025
    Dil_const::Float64 = 0.072
    threshold::Float64 = 10.0
end

@with_kw struct ini_conditions
    N₁⁰::Float64 = 10.0
    N₂⁰::Float64 = 0.2
    C⁰::Float64 = 0.01
    SC₁⁰::Float64 = 0.1
    SC₂⁰::Float64 = 0.1
    Δ⁰::Int = 0
end

# Set the default values of the parameters and the initial conditions
model_constants = parameters()
u0 = ini_conditions()

# make functions one liners for readability
sub_asm(a, b, X) = a * X / (b + X)
dilution(X, Y, Z, a, b, c, d) = a * X + b * Y + c * Z + d

# Define the ODE system, parameter list p can be retrieved from the previous settings
function ode_sys_tri(du, u, p, t)
    du[1] =
        u[6] * dilution(u[3], u[4], u[5], p.pop_c, p.pop_sc_1, p.pop_sc_2, 0) * (p.N₁⁰ - u[1]) -
        p.c₁ * sub_asm(p.a, p.b, u[1]) * u[3] / p.OD_chlamy -
        p.c₂ * sub_asm(p.c_1, p.d_1, u[1]) * u[4] / p.OD_bacteria - 
        p.c₄ * sub_asm(p.c_2, p.d_2, u[1]) * u[5] / p.OD_bacteria

    du[2] =
        p.α * u[3] -
        p.c₃ * sub_asm(p.e_1, p.f_1, u[2]) * u[4] / p.OD_bacteria -
        p.c₅ * sub_asm(p.e_2, p.f_2, u[2]) * u[5] / p.OD_bacteria -
        u[6] * u[2] * dilution(u[3], u[4], u[5], p.pop_c, p.pop_sc_1, p.pop_sc_2, 0)

    du[3] =
        u[3] * (
            p.c₁ * sub_asm(p.a, p.b, u[1]) -
            u[6] * dilution(u[3], u[4], u[5], p.pop_c, p.pop_sc_1, p.pop_sc_2, 0) - 
            p.d₁
        )

    du[4] =
        u[4] * (
            p.c₂ * sub_asm(p.c_1, p.d_1, u[1]) +
            p.c₃ * sub_asm(p.e_1, p.f_1, u[2]) - p.d₂ -
            u[6] * dilution(u[3], u[4], u[5], p.pop_c, p.pop_sc_1, p.pop_sc_2, 0)
        )

    du[5] = 
          u[5] * (
            p.c₄ * sub_asm(p.c_2, p.d_2, u[1]) +
            p.c₅ * sub_asm(p.e_2, p.f_2, u[2]) - p.d₃ -
            u[6] * dilution(u[3], u[4], u[5], p.pop_c, p.pop_sc_1, p.pop_sc_2, 0)
        )

    du[6] = 0
end

# Modify the conditions for the dilution switch
activate_dilution_condition(u, t, integrator) =
    integrator.u[6] < 1 && integrator.u[3] + integrator.u[4] + integrator.u[5] > integrator.p.threshold

activate_dilution_affect!(integrator) = integrator.u[6] = 1

stop_dilution_condition(u, t, integrator) =
    integrator.u[6] > 0 && integrator.u[3] + integrator.u[4] + integrator.u[5] < 0.9 * integrator.p.threshold

stop_dilution_affect!(integrator) = integrator.u[6] = 0.0

# Create the callback set with switching conditions
cb1 = DiscreteCallback(activate_dilution_condition, activate_dilution_affect!)
cb2 = DiscreteCallback(stop_dilution_condition, stop_dilution_affect!)
cbs = CallbackSet(cb1, cb2)

# Set the time for simulation
tspan = (0.0, 1000.0)

# Define the ODE problem for simulation
prob = ODEProblem(ode_sys, u0, tspan, model_constants)
sol = solve(
    prob,
    callback=cbs,
    saveat=0.1
)

plot_layout = @layout [
    grid(3,2)
    a{0.2h}
]

plot(sol,
    dpi = 500,
    size = (800,600),
    layout = plot_layout,
    title = ["N₁" "N₂" "Chlamy" "SynCom_1" "SynCom_2" "Dilution Switch"],
    label = ["" "" "" "" "" ""]
    )

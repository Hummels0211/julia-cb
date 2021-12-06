using DifferentialEquations
using Plots

# Set the initial conditions
N₁⁰ = 10.0
N₂⁰ = 0.2
C⁰ =  0.01
SC⁰ = 0.1

# Parameter table
#############################################################
c₁ = 0.1      # constant for substrate intake by chlamy
c₂ = 0.35      # constant for substrate intake by bacteria
c₃ = 0.31     # constant for vitamin intake by bacteria
#############################################################
a = 1.5
b = 7.5        # constants for chlamy intake function
c = 0.34
d = 10.2       # constants for bacteria intake function
e = 0.37
f = 0.52        # constants for bacteria vitamin intake function
#############################################################
# The population effect for two species (optional)
OD_chlamy = 1
OD_bacteria = 1
#############################################################
α = 2     # constant for carbon production by chlamy
g₁ = 10      # constant for chlamy growth limiting effect
g₂ = 12      # constant for bacteria growth limiting effect

d₁ = 0.005      # constant for chlamy death
d₂ = 0.05      # constant for bacterial death
#############################################################
pop_c = 0.06
pop_sc = 0.06    # population constant
#############################################################

# Setting the threshold of the biomass
threshold = 10

# Initialising the dilution states to 0
Δ = 0

# Setting the dilution rate of the media
Dil_const = 0.072

################################################################################
# Define the function for substrate intake by the microbial organisms
################################################################################
# The default Monod function
function sub_asm(a, b, X)
    asm = a * X / (b + X)
    return asm
end

# The linear assimilation pattern
function sub_asm_lin(a, X)
    asm = a * X
    return asm
end

# The Michaelis-Menten assimilation pattern
function sub_sam_cat(a, b, X)
    asm = a * (X^2) / (b^2 + X^2)
    return asm
end

# The exponential assimilation pattern
function sub_sam_exp(a, b, X)
    asm = a * (1 - np.exp(- b * X))
    return asm
end

################################################################################

# Define the growth function
function growth_sat(a, X)
    growth = a/(a + X)
    return growth
end

# Define the dilution function under turbidostat condition
function dilution(X, Y, a, b, d)
    return a * X + b * Y + d
end


# Define the ODE system, parameter list p can be retrieved from the previous settings
function ode_sys(du, u, p, t)
    du[1] = Δ * dilution(u[3], u[4], pop_c, pop_sc, 0) * (N₁⁰ - u[1]) -
            c₁ * sub_asm(a, b, u[1]) *u[3] / OD_chlamy  -
            c₂ * sub_asm(c, d, u[1]) * u[3] / OD_bacteria
    du[2] = α * u[3] - c₃ * sub_asm(e, f, u[2]) * u[4] / OD_bacteria -
            Δ * u[2] * dilution(u[3], u[4], pop_c, pop_sc, 0)
    du[3] = u[3] * (c₁ * sub_asm(a, b, u[1]) -
            Δ * dilution(u[3], u[4], pop_c, pop_sc, 0) - d₁)
    du[4] = u[4] + u[4] * (c₂ * sub_asm(c, d, u[1]) +
            c₃ * sub_asm(e, f, u[2] - d₂) - Δ * dilution(u[3], u[4], pop_c, pop_sc, 0))
end

# Set the condition of trigging dilution in the bioreactor
function condition_dilute(u, t, integrator)
    # The condition when the biomass of chlamydomonas & bacteria surpass the threshold
    integrator.u[3] + integrator.u[4] > threshold
end

function affect_dilute!(integrator)
    Δ = 1  # Turn on the dilution
end

# Set the condition of stopping dilution in the bioreactor
function condition_off(u, t, integrator)
    # The dilution should be stopped when the biomass reaches the lower boundary
    integrator.u[3] + integrator.u[4] <= 0.9 * threshold
end

function affect_stop!(integrator)
    Δ = 0 # Stop diluting the bioreactor, restore to closed system
end

# Set the call back sets
dilut_cb = ContinuousCallback(condition_dilute, affect_dilute!)
stop_cb = ContinuousCallback(condition_off, affect_stop!)
cb = CallbackSet(dilut_cb, stop_cb)

# Set the timespan for simulation
tspan = (0.0, 20.0)


# Set the initial condition
u0 = [N₁⁰, N₂⁰, C⁰, SC⁰]

# Solve the ODE system (with different solvers [optional])
prob = ODEProblem(ode_sys, u0, tspan)
sol = solve(prob, DP5(), callback = cb)


plot(sol)

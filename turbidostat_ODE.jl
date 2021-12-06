using DifferentialEquations
using Plots

# Define the constants
c_1 = 1      # constant for substrate intake by chlamy
c_2 = 1      # constant for substrate intake by bacteria
C = 2        # constant for carbon production by chlamy
g_1 = 0.1      # constant for chlamy growth limiting effect
g_2 = 1      # constant for bacteria growth limiting effect
d_1 = 0.001      # constant for chlamy death
d_2 = 0.001      # constant for bacterial death
a = 0.25
b = 0.5        # constants for chlamy intake function
c = 0.25
d = 0.25       # constants for bacteria intake function

# Setting the dilution rate of the media
Dil_const = 0.2


# Define the function for substrate intake by the microbial organisms
function sub_asm(a, b, X)
    asm = 1 - a/(X + b)
    return asm
end

# Define the growth function
function growth_sat(a, X)
    growth = a / (a + X)
    return growth
end

# Define the
function simp_model(d_var, var, t, p)
    c_1, c_2, C, g_1, g_2, d_1, d_2, a, b, c, d = p
    if var[3] + var[4] <= threshold
        D = 0
    else
        D = Dil_const
    end
    d_var[1] = - c_1 * sub_asm(a, b, var[3]) + D
    d_var[2] = - c_2 * sub_asm(c, d, var[4]) + C * var[3] - D
    d_var[3] = growth_sat(g_1, var[3]) * sub_asm(a, b, var[3]) - (d_1 + D) * var[3]
    d_var[4] = growth_sat(g_2, var[4]) * sub_asm(c, d, var[4]) - (d_2 + D) * var[4]
end


tspan = (0.0,72.0)
var_0 = [10.0, 0.0, 0.02, 0.0002]
prob = ODEProblem(simp_model, tspan, var_0)

sol = solve(prob)


plot(sol)

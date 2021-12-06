using DifferentialEquations
using Plots

a = 0.02
b = 0.2
c = -65
d = 8
i = 0

p = [a,b,c,d,i]

function fun(du,u,p,t)
    if u[1] <30
        du[1] = (0.04*u[1] + 5)*u[1] + 150 - u[2] - p[5]
        du[2] = p[1]*(p[2]*u[1]-u[2])
    else
        u[1] = p[3]
        u[2] = u[2] + p[4]
    end
end

u0 = [0.0;0.0]
tspan = (0.0,100)
prob = ODEProblem(fun,u0,tspan,p)
# tic()
sol = solve(prob,reltol = 1e-8)
# toc()

plot(sol)

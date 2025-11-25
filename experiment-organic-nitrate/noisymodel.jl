using DifferentialEquations
using LinearAlgebra
using Plots
using Random
using Statistics

# seed the random number generator for reproducibility
Random.seed!(1)

# coupled ODE system of chemical species

#                     Photolysis rate constant	Hydrolysis rate constant	Alkoxy radical branching ratio	Nitrate branching ratio from RO2 + NO
# Phase	              j1 (s–1)	                k8 (s–1)	                k3[O2]/(k3[O2] + k2)	        X5a
# Gas	              (2.9 ± 0.1) × 10–4	    N/A	                        (6.6 ± 0.4) × 10–3	            0.22 ± 0.001
# Organic	          (4.3 ± 0.6) × 10–3	    N/A	                        0.099 ± 0.008	                0.86 ± 0.04
# Aqueous	          (7.0 ± 0.2) × 10–5	    (1.6 ± 0.7) × 10–3	        0.31 ± 0.007	                0.79 ± 0.03

function barber2024!(du, u, p, t)
    # Unpack species concentrations
    A, B, D, G, IplusJ = u
    # Unpack rate constants
    j1, k8, R, X5a = p
    
    # Define the ODEs
    du[1] = -j1 * A  - k8 * A    # dA/dt  /\/\/ONO, n-pentyl nitrite precursor colored brown
    du[2] = k8 * A               # dB/dt  C5H12O, colored green
    du[3] = R/(R+1)*j1*A         # dD/dt  C5H10O, colored goldenrod
    du[4] = X5a*(1/(R+1))*j1*A   # dG/dt  C5H11NO4, colored blue 
    du[5] = (1-X5a)/(R+1)*j1*A   # d([I]+[J])/dt C5H10O2, colored red
end

alkoxy_radical_branching_ratios = [6.6e-3, 0.099, 0.31]
# note R is defined as k3[O2]/k2, or kRO+O2/kH-shift in the supplement
# so B = R/(1+R) -> R = B/(1-B)
R = (alkoxy_radical_branching_ratios)./(1 .- alkoxy_radical_branching_ratios)


p_gas =     [(2.9e-4), 0.0,     R[1], 0.22]
p_organic = [(4.3e-3), 0.0,     R[2], 0.86]
p_aqueous = [(7.0e-5),(1.6e-3), R[3], 0.79]

# Initial concentrations
A0 = 45           # Initial concentration of A
B0 = 0.0          # Initial concentration of B
D0 = 0.0          # Initial concentration of D
G0 = 0.0          # Initial concentration of G
IplusJ0 = 0.0     # Initial concentration of [I] + [J]
u0 = [A0, B0, D0, G0, IplusJ0]

dt = 1.0
tspan = (0.0, 1200.0)  # Time span for the simulation

# specify solver
# solver = Rosenbrock23()
solver = QNDF() # analogous to MATLAB's ode15s used by Tori

# Solve for gas phase
prob_gas = ODEProblem(barber2024!, u0, tspan, p_gas)
sol_gas = solve(prob_gas, solver, saveat=dt)
# Solve for organic phase
prob_organic = ODEProblem(barber2024!, u0, tspan, p_organic)
sol_organic = solve(prob_organic, solver, saveat=dt)
# Solve for aqueous phase
prob_aqueous = ODEProblem(barber2024!, u0, tspan, p_aqueous)
sol_aqueous = solve(prob_aqueous, solver, saveat=dt)

# plot a stacked subplot for the three experiments
# colors: A - brown, B - green, D - goldenrod, G - blue, I+J - red
using Plots
p = plot(layout = (3, 1), size=(500, 700), grid=false)


# gas phase on subplot 1
plot!(p[1], sol_gas.t, sol_gas[1, :], label="PN precursor", color=:brown)
plot!(p[1], sol_gas.t, sol_gas[2, :], label="\$C_5H_{12}O\$", color=:green)
plot!(p[1], sol_gas.t, sol_gas[3, :], label="\$C_5H_{10}O\$", color=:goldenrod)
plot!(p[1], sol_gas.t, sol_gas[4, :], label="\$C_5H_{11}NO_4\$", color=:blue)
plot!(p[1], sol_gas.t, sol_gas[5, :], label="\$C_5H_{10}O_2\$", color=:red)
title!(p[1], "Gas Phase Experiment")
xlabel!(p[1], "Time (s)")
ylabel!(p[1], "Concentration")

# organic phase on subplot 2
plot!(p[2], sol_organic.t, sol_organic[1, :], label="PN precursor", color=:brown)
plot!(p[2], sol_organic.t, sol_organic[2, :], label="\$C_5H_{12}O\$", color=:green)
plot!(p[2], sol_organic.t, sol_organic[3, :], label="\$C_5H_{10}O\$", color=:goldenrod)
plot!(p[2], sol_organic.t, sol_organic[4, :], label="\$C_5H_{11}NO_4\$", color=:blue)
plot!(p[2], sol_organic.t, sol_organic[5, :], label="\$C_5H_{10}O_2\$", color=:red)
title!(p[2], "Organic Phase Experiment")
xlabel!(p[2], "Time (s)")
ylabel!(p[2], "Concentration")  

# aqueous phase on subplot 3
plot!(p[3], sol_aqueous.t, sol_aqueous[1, :], label="PN precursor", color=:brown)
plot!(p[3], sol_aqueous.t, sol_aqueous[2, :], label="\$C_5H_{12}O\$", color=:green)
plot!(p[3], sol_aqueous.t, sol_aqueous[3, :], label="\$C_5H_{10}O\$", color=:goldenrod)
plot!(p[3], sol_aqueous.t, sol_aqueous[4, :], label="\$C_5H_{11}NO_4\$", color=:blue)
plot!(p[3], sol_aqueous.t, sol_aqueous[5, :], label="\$C_5H_{10}O_2\$", color=:red)
title!(p[3], "Aqueous Phase Experiment")
xlabel!(p[3], "Time (s)")
ylabel!(p[3], "Concentration")


data_gas = hcat(sol_gas[1, :], sol_gas[2, :], sol_gas[3, :], sol_gas[4, :], sol_gas[5, :])
data_organic = hcat(sol_organic[1, :], sol_organic[2, :], sol_organic[3, :], sol_organic[4, :], sol_organic[5, :])
data_aqueous = hcat(sol_aqueous[1, :], sol_aqueous[2, :], sol_aqueous[3, :], sol_aqueous[4, :], sol_aqueous[5, :])

# add % noise to each solution for PCA
percent_noise = 0.10
noise_gas = 1 .+ percent_noise .* randn(size(data_gas))
noise_organic = 1 .+ percent_noise .* randn(size(data_organic))
noise_aqueous = 1 .+ percent_noise .* randn(size(data_aqueous))

noisy_data_gas = data_gas .* noise_gas
noisy_data_organic = data_organic .* noise_organic
noisy_data_aqueous = data_aqueous .* noise_aqueous

# smoothing function
function smooth_data(data, window_size)
    smoothed_data = similar(data)
    half_window = div(window_size, 2)
    for i in 1:size(data, 2)
        for j in 1:size(data, 1)
            start_idx = max(1, j - half_window)
            end_idx = min(size(data, 1), j + half_window)
            smoothed_data[j, i] = mean(data[start_idx:end_idx, i])
        end
    end
    return smoothed_data
end
# smooth noisy data with a moving average filter of window size 15
window_size = 15
noisy_data_gas = smooth_data(noisy_data_gas, window_size)
noisy_data_organic = smooth_data(noisy_data_organic, window_size)
noisy_data_aqueous = smooth_data(noisy_data_aqueous, window_size)

# add noisy points to the first plot p[1], p[2], p[3]
plot!(p[1], sol_gas.t, noisy_data_gas[:,1], seriestype=:scatter,label=false, color=:brown, markerstrokecolor=:brown, markersize=0.5, alpha=0.8)
plot!(p[1], sol_gas.t, noisy_data_gas[:,2], seriestype=:scatter, label=false, color=:green, markerstrokecolor=:green, markersize=0.5, alpha=0.8)
plot!(p[1], sol_gas.t, noisy_data_gas[:,3], seriestype=:scatter, label=false, color=:goldenrod, markerstrokecolor=:goldenrod, markersize=0.5, alpha=0.8)
plot!(p[1], sol_gas.t, noisy_data_gas[:,4], seriestype=:scatter, label=false, color=:blue, markerstrokecolor=:blue, markersize=0.5, alpha=0.8)
plot!(p[1], sol_gas.t, noisy_data_gas[:,5], seriestype=:scatter, label=false, color=:red, markerstrokecolor=:red, markersize=0.5, alpha=0.8)
plot!(p[2], sol_organic.t, noisy_data_organic[:,1], seriestype=:scatter,label=false, color=:brown, markerstrokecolor=:brown, markersize=0.5, alpha=0.8)
plot!(p[2], sol_organic.t, noisy_data_organic[:,2], seriestype=:scatter, label=false, color=:green, markerstrokecolor=:green, markersize=0.5, alpha=0.8)
plot!(p[2], sol_organic.t, noisy_data_organic[:,3], seriestype=:scatter, label=false, color=:goldenrod, markerstrokecolor=:goldenrod, markersize=0.5, alpha=0.8)
plot!(p[2], sol_organic.t, noisy_data_organic[:,4], seriestype=:scatter, label=false, color=:blue, markerstrokecolor=:blue, markersize=0.5, alpha=0.8)
plot!(p[2], sol_organic.t, noisy_data_organic[:,5], seriestype=:scatter, label=false, color=:red, markerstrokecolor=:red, markersize=0.5, alpha=0.8)
plot!(p[3], sol_aqueous.t, noisy_data_aqueous[:,1], seriestype=:scatter,label=false, color=:brown, markerstrokecolor=:brown, markersize=0.5, alpha=0.8)
plot!(p[3], sol_aqueous.t, noisy_data_aqueous[:,2], seriestype=:scatter, label=false, color=:green, markerstrokecolor=:green, markersize=0.5, alpha=0.8)
plot!(p[3], sol_aqueous.t, noisy_data_aqueous[:,3], seriestype=:scatter, label=false, color=:goldenrod, markerstrokecolor=:goldenrod, markersize=0.5, alpha=0.8)
plot!(p[3], sol_aqueous.t, noisy_data_aqueous[:,4], seriestype=:scatter, label=false, color=:blue, markerstrokecolor=:blue, markersize=0.5, alpha=0.8)
plot!(p[3], sol_aqueous.t, noisy_data_aqueous[:,5], seriestype=:scatter, label=false, color=:red, markerstrokecolor=:red, markersize=0.5, alpha=0.8)
display(p)

# save figure as high resolution png
savefig(p, "./experiment-organic-nitrate/noisymodel.pdf")

# # run PCA on each solution, plot all 5 principal components for each phase in a single plot
using MultivariateStats
Mgas = fit(PCA, noisy_data_gas'; pratio=1.0-1e-24)
Morganic = fit(PCA, noisy_data_organic'; pratio=1.0-1e-24)
Maqueous = fit(PCA, noisy_data_aqueous'; pratio=1.0-1e-24)

# plot the first 5 eigvals for each phase each as a line
p2 = plot(size=(600,400),grid=false)
plot!(p2, eigvals(Mgas)[1:4], label="Gas Phase", marker=:circle, yscale=:log10)
plot!(p2, eigvals(Morganic)[1:4], label="Organic Phase", marker=:square, yscale=:log10)
plot!(p2, eigvals(Maqueous)[1:4], label="Aqueous Phase", marker=:diamond, yscale=:log10)
xlabel!(p2, "Principal Component")
ylabel!(p2, "Eigenvalue")
display(p2)

# save figure as png 
# png(p2, "figure2_eigenvalues.png")


eigvecs(Mgas)[:,1]/(-eigvecs(Mgas)[1,1])
eigvecs(Morganic)[:,1]/(-eigvecs(Morganic)[1,1])
eigvecs(Maqueous)[:,1]/(-eigvecs(Maqueous)[1,1])

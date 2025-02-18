using DifferentialEquations, SpecialFunctions, LinearAlgebra, Plots, QuadGK,SphericalFunctions

lmax = 10; 
l_num = (lmax +1)^2

# Define a small range of θ values
theta = range(0, π, length=100);
phi = range(0, 2π, length=100);

# Define an initial condition in terms of spherical harmonics
function gauss(theta)
    return exp(-(theta - 0.01)^2/(2*0.2^2))  # A Gaussian decay in l-space
end

for l in 0:lmax
    for m in -l:l
        if m == 0  # Only consider the m = 0 mode
          # Define spherical harmonic function Y_lm for m = 0, varying l
            Ylm = theta -> real(sYlm_values(theta, 0.0, l, m)[1])  # Real part of Y_lm

            # Define the integrand function
            integrand(theta) = gauss(theta) * Ylm(theta) * sin(theta)

            # Compute the coefficient c_l using numerical integration (quadgk)
            c_l, _ = quadgk(integrand, 0, π)
            # Append the computed coefficient to the list
            push!(c_i, c_l)
        else
            push!(c_i,0)
        end
    end
end

# function for the A matrix
function A_!(dA, A, p, t)
    L_ = p[1]  # Laplacian eigenvalues 
    
    # Time derivatives
    n = length(A) ÷ 2 
    clm = A[1:n]       
    dt_c = A[n+1:end] 

    @. dA[1:n] = dt_c 
    @. dA[n+1:end] = L_ * clm  
end

# Arrays for the 
c_i = zeros(Float64, l_num)
dtc_i = zeros(Float64, l_num)

# Example 
c_i[2] = 1.0 
u_i = vcat(c_i, dtc_i) 

# Define time span for ODE solver
tspan = (0.0, 20.0)  

# Define the ODE problem, passing Lp_eigen as a parameter
prob = ODEProblem(A_!, u_i, tspan, [Lp_eigen])

# 4th order Runge Kutta solver
sol = solve(prob, Tsit5(), reltol=1e-6, abstol=1e-6)

# Plotting: Mode (l=1, m=0)
plot(sol.t, sol[2, :], xlabel="Time", ylabel="Amplitude", title="Mode (l=1, m=0)")



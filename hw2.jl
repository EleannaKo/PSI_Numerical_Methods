using DifferentialEquations, SpecialFunctions, LinearAlgebra, Plots, QuadGK, SphericalFunctions

# Setting up the problem: Wave equation on a sphere

# Define the maximum spherical harmonic degree
lmax = 10
ltot = (lmax + 1)^2  # Total spherical harmonic modes

# Compute eigenvalues of the Laplacian operator
eigenvalues = Float64[]
for l in 0:lmax
    for m in -l:l
        push!(eigenvalues, -l * (l + 1))
    end
end

# A_matrix function
function A_matrix!(dA, A, params, t)
    # Extract Laplacian eigenvalues
    eigenvalues = params[1]
    # Determine number of coefficients
    half_n = length(A) ÷ 2
    coeffs = A[1:half_n]
    velocities = A[half_n+1:end]

    # Compute derivatives
    @. dA[1:half_n] = velocities
    @. dA[half_n+1:end] = eigenvalues * coeffs
end

# Initialize arrays for coefficients and their time derivatives
coeff_arr = zeros(Float64, ltot)
velocity_arr = zeros(Float64, ltot)

# Setting initial condition: Exciting (l=1, m=0) mode
coeff_arr[2] = 1.0  # Nonzero initial value for (l=1, m=0)
init_state = vcat(coeff_arr, velocity_arr)  # Merge arrays

# Define the simulation time range
time_range = (0.0, 20.0)

# Set up and solve the ODE system
problem = ODEProblem(A_matrix!, init_state, time_range, [eigenvalues])
solution = solve(problem, Tsit5(), reltol=1e-6, abstol=1e-6)

# Plot the results for (l=1, m=0) mode
plot(solution.t, solution[2, :], xlabel="T", ylabel="A", title="l=1, m=0")

# Gaussian function: Computes the peak value of a Gaussian curve
std_dev = 0.2  
function gaussian_peak(angle)
    return exp(-((angle - 0.001)^2) / (2 * std_dev^2))
end

# Computing mode coefficients via numerical integration
function compute_coefficients(lmax)
    coeff_values = Float64[]  # Array to store calculated coefficients
    
    for l in 0:lmax
        for m in -l:l
            if m == 0  # Focusing only on m = 0
                Y_lm_func = θ -> real(sYlm_values(θ, 0.0, l, m)[1])  # Extract real part
                
                # Define the function to be integrated
                integral_func(θ) = gaussian_peak(θ) * Y_lm_func(θ) * sin(θ)
                
                # Perform numerical integration
                coeff, _ = quadgk(integral_func, 0, π)
                push!(coeff_values, coeff)
            else
                push!(coeff_values, 0)
            end
        end
    end
    
    return coeff_values  # Return computed coefficients
end

# Example setup: Solving wave equation
max_degree = 10
num_spherical_modes = (max_degree + 1)^2  # Compute total spherical harmonic modes

# Obtain initial coefficient values
initial_coeffs = compute_coefficients(max_degree)

# Initialize time derivatives (zero except for a specific mode)
deriv_coeffs = zeros(Float64, num_spherical_modes)
deriv_coeffs[4] = 0.1  # Introduce excitation in (l=1, m=0)

# Merge coefficients and derivatives into a single vector
initial_state = vcat(initial_coeffs, deriv_coeffs)

# Define time range for solving the system
simulation_time = (0.0, 10.0)

# Construct and solve the differential equation model
wave_problem = ODEProblem(A_matrix!, initial_state, simulation_time, [eigenvalues])
wavesolution = solve(wave_problem, Tsit5(), reltol=1e-6, abstol=1e-6)

# Visualization: Plot evolution of the (l=1, m=0) mode
plot(wavesolution.t, wavesolution[4, :], xlabel="T", ylabel="A", 
     title="Gaussian initial condition peaking at the North Pole")





# Initialize theta, phi angles
theta_ = range(0, π, length=100)
phi_ = range(0, 2π, length=100)

# Create the meshgrid
theta = repeat(theta_, 1, length(phi_))  # Repeat theta along columns
phi = repeat(phi_', length(theta_), 1)   # Repeat phi along rows

# Compute Cartesian coordinates using broadcasting
x = sin.(theta) .* cos.(phi)
y = sin.(theta) .* sin.(phi)
z = cos.(theta)

# Convert into 1D arrays for scatter plotting
x_flat = vec(x)
y_flat = vec(y)
z_flat = vec(z)

# Plot the spherical mesh
scatter3d(x_flat, y_flat, z_flat, marker = (:auto, 3), 
          xlabel = "x", ylabel = "y", zlabel = "z", 
          title = "Spherical Grid")





# Reconstruct the scalar field from solution coefficients
function r_psi(solution, time_idx, lmax, theta_grid, phi_grid)
    psi_field = zeros(size(theta_grid))  # Initialize the field
    
    coeff_idx = 1  # Index for accessing solution coefficients
    
    for l in 0:lmax
        for m in -l:l
            # Compute real spherical harmonics over the grid
            Ylm_grid = [real(sYlm_values(theta_grid[row, col], phi_grid[row, col], l, m)[1]) 
                        for row in 1:size(theta_grid, 1), col in 1:size(theta_grid, 2)]
            
            # Retrieve the coefficient from the solution array
            coeff_lm = solution[coeff_idx, time_idx]
            
            # Accumulate contributions from different modes
            @. psi_field += coeff_lm * Ylm_grid  
            
            coeff_idx += 1  # Move to the next coefficient
        end
    end
    
    return psi_field
end

animation = @animate for time_idx in 1:30:length(solution.t)
    psi_field = r_psi(solution, time_idx, 2, theta, phi)
    
    heatmap(phi_, theta_, psi_field, 
            title="Wave Evolution at t = $(round(solution.t[time_idx], digits=2))",
            ylabel="ϕ", xlabel="θ", clims=(-0.1, 0.1))
end
    
# Export the animation as a GIF
gif(animation, "lmax2.gif", fps=5)


animation = @animate for time_idx in 1:30:length(solution.t)
    psi_field = r_psi(solution, time_idx, 5, theta, phi)
    
    heatmap(phi_, theta_, psi_field,
            title="Wave Evolution at t = $(round(solution.t[time_idx], digits=2))",
            ylabel="ϕ", xlabel="θ", clims=(-0.1, 0.1))
end

# Export the animation as a GIF
gif(animation, "lmax5.gif", fps=5)


animation = @animate for time_idx in 1:30:length(solution.t)
    psi_field = r_psi(solution, time_idx, 10, theta, phi)
    
    heatmap(phi_, theta_, psi_field,  
            title="Wave Evolution at t = $(round(solution.t[time_idx], digits=2))",
            ylabel="ϕ", xlabel="θ", clims=(-0.1, 0.1))
end

# Export the animation as a GIF
gif(animation, "lmax10.gif", fps=5)
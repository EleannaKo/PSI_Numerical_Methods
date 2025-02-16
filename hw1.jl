using Pkg
Pkg.add("Interpolations")
using WGLMakie
using Interpolations


# Let's create the data set of xx,yy -- I'll use y = sin(x) as an example.3
min_x = 0.01
max_x = 6.28
xx = range(min_x, max_x, 10)
yy = sin.(xx);
# This is the weird syntax of this "Interpolations" package for a cubic spline, on a grid, with natural boundaries...
spl = interpolate(yy, BSpline(Cubic(Natural(OnGrid()))))
spl = scale(spl, xx);

# f = Figure()
# Axis(f[1,1], title="Data and Spline", xlabel="X", ylabel="Y")
# scatter!(xx, yy)
# plotx = range(min_x, max_x, 1000)
# lines!(plotx, spl(plotx), linewidth=3)
# lines!(plotx, sin.(plotx))
# Axis(f[2,1], title="Error in spline", xlabel="X", ylabel="delta-Y")
# lines!(plotx, spl(plotx) .- sin.(plotx))
# scatter!(xx, zeros(length(xx)))
# f

# Our data points go from x_0 to x_N, so N is actually one less than the number of data points in xx!!

# With Julia's 1-based indexing, things are going to be a bit messy, because x_0 is actually in xx[1] !!

N = length(xx) - 1
A = zeros((4*N, 4*N))
r = zeros(4*N);

# .... YOUR CODE HERE to fill in elements of the A matrix

# This is tricky... it took me more than a few tries to get it right!

# You'll definitely want to write out the equations for each of the constraints (on paper)
# before trying to write the code here.

# Note that the coefficient will be laid out so that, in the matrix A,
#   terms that multiply a_i will be in column 4*(i-1) + 1
#                       b_i                   4*(i-1) + 2
#                       c_i                   4*(i-1) + 3
#                       d_i                   4*(i-1) + 4
# and if you need to refer to the coefficients for i+1, those will be in
#                       a_{i+1}               4*(i-1) + 5
#                       b_{i+1}               4*(i-1) + 6
#                       c_{i+1}               4*(i-1) + 7
#                       d_{i+1}               4*(i-1) + 8

# Left endpoint constraints: C_i(x_{i-1}) = y_{i-1} or in julia notation C_i(x_{i}) = y_{i}
# Put these in the first N rows...
for i in 1:N  # N here corresponds to the polynomial C[i] where i goes from 1 to N 
    A[i, (4*(i-1))+1] = 1
    A[i, (4*(i-1))+2] = xx[i]
    A[i, (4*(i-1))+3] = xx[i]^2
    A[i, (4*(i-1))+4] = xx[i]^3
    r[i] = yy[i]
end
    
    # Right endpoint constraints: C_i(x_i) = y_i or in julia notation C_i(x_{i+1}) = y_{i+1}
    # Put these in the next N rows, N+1 to 2*N
for i in 1:N
    A[N + i, (4*(i-1))+1] = 1
    A[N + i, (4*(i-1))+2] = xx[i+1]
    A[N + i, (4*(i-1))+3] = xx[i+1]^2
    A[N + i, (4*(i-1))+4] = xx[i+1]^3
    r[N + i] = yy[i+1]
end

    # First derivative constraints: C_i'(x_i) = C_{i+1}'(x_i) or in julia notation C_i'(x_{i+1}) = C_{i+1}'(x_{i+1})
    # (note that this only goes up to N-1)
for i in 1:N-1
    A[2*N + i, (4*(i-1))+1] = 0
    A[2*N + i, (4*(i-1))+2] = 1
    A[2*N + i, (4*(i-1))+3] = 2*xx[i+1]
    A[2*N + i, (4*(i-1))+4] = 3*xx[i+1]^2
    A[2*N + i, (4*(i-1))+5] =-A[2*N + i, (4*(i-1))+1]
    A[2*N + i, (4*(i-1))+6] =-A[2*N + i, (4*(i-1))+2]
    A[2*N + i, (4*(i-1))+7] =-A[2*N + i, (4*(i-1))+3]
    A[2*N + i, (4*(i-1))+8] =-A[2*N + i, (4*(i-1))+4]
    r[2*N + i] = 0
end
    
    # Second derivative constraints: C_i''(x_i) = C_{i+1}''(x_i) or in julia notation C_i"(x_{i+1}) = C_{i+1}"(x_{i+1})
for i in 1:N-1
    A[3*N + i, (4*(i-1))+1] = 0
    A[3*N + i, (4*(i-1))+2] = 0
    A[3*N + i, (4*(i-1))+3] = 2
    A[3*N + i, (4*(i-1))+4] = 6*xx[i+1]
    A[3*N + i, (4*(i-1))+5] =-A[3*N + i, (4*(i-1))+1]
    A[3*N + i, (4*(i-1))+6] =-A[3*N + i, (4*(i-1))+2]
    A[3*N + i, (4*(i-1))+7] =-A[3*N + i, (4*(i-1))+3]
    A[3*N + i, (4*(i-1))+8] =-A[3*N + i, (4*(i-1))+4]
    r[3*N + i] = 0
end
    
    # Natural constraints for the leftmost and rightmost data point...
    # You can put these in rows 3*N and 4*N, because those aren't filled in yet!
 A[3*N, 3] = 2;
 A[3*N, 4] = 6*xx[1];
 
 A[4*N, (4*(N-1))+3] = 2;
 A[4*N, (4*(N-1))+4] = 6*xx[N+1];

 # Now we get our coefficient vector by solving the matrix equation:
c = A \ r;

# One useful feature for debugging is the spy() function -- this shows you which elements of a matrix are
# non-zero.
# f = Figure()
# Axis(f[1,1], title="A matrix elements", ylabel="ROW", xlabel="COLUMN")
# spy!(A')
# f
     

# # For debugging, you might find it helpful to plot the first few function C_1, C_2, C_3. 

# f = Figure()
# ax = Axis(f[1,1], title="Data and Spline", xlabel="X", ylabel="Y")
# scatter!(xx, yy)
# plotx = range(min_x, max_x, 1000)

# C1 = x -> @. c[1] + c[2]  * x + c[3]  * x^2 + c[4]  * x^3
# C2 = x -> @. c[5] + c[6]  * x + c[7]  * x^2 + c[8]  * x^3
# C3 = x -> @. c[9] + c[10] * x + c[11] * x^2 + c[12] * x^3
# margin = 0
# lines!(range(xx[1]-margin, xx[2]+margin, 20), C1, label="C1")
# lines!(range(xx[2]-margin, xx[3]+margin, 20), C2, label="C2")
# lines!(range(xx[3]-margin, xx[4]+margin, 20), C3, label="C3")
# ylims!(-2, +2)
# f[1,2] = Legend(f, ax)
# f

# Here's one way we can create a single function to evaluate the spline:
# This is a function that creates and returns another function!
function make_spline(xgrid, coeffs)
    # Here's the spline function we're creating...
    function spline_function(x)
        N = length(xgrid)-1
        for i in 1:N
            # Is this the right segment?
            if (x < xgrid[i+1]) || (i == N)
                # Evaluate the cubic polynomial for this segment!
                i0 = 4*(i-1)
                return coeffs[i0 + 1] + coeffs[i0 + 2] * x + coeffs[i0 + 3] * x^2 + coeffs[i0 + 4] * x^3
            end
        end
    end
    # We can use the "map" function so that a scalar input produces a scalar output,
    # and a vector input produces a vector output
    function vectorized_spline_function(x)
        return map(spline_function, x)
    end
    return vectorized_spline_function
end;

spline = make_spline(xx, c);
spline(5.)
spline([5., 6.])
f = Figure()
ax = Axis(f[1,1], title="Data and Spline", xlabel="X", ylabel="Y")
scatter!(xx, yy)
plotx = range(min_x, max_x, 1000)
lines!(plotx, spline)
f

################################################################

function get_spline_coefficients(xx, yy)
    N = length(xx) - 1
    A = zeros((4*N, 4*N))
    r = zeros(4*N)
    for i in 1:N  # N here corresponds to the polynomial C[i] where i goes from 1 to N 
        A[i, (4*(i-1))+1] = 1
        A[i, (4*(i-1))+2] = xx[i]
        A[i, (4*(i-1))+3] = xx[i]^2
        A[i, (4*(i-1))+4] = xx[i]^3
        r[i] = yy[i]
    end
        
        # Right endpoint constraints: C_i(x_i) = y_i or in julia notation C_i(x_{i+1}) = y_{i+1}
        # Put these in the next N rows, N+1 to 2*N
    for i in 1:N
        A[N + i, (4*(i-1))+1] = 1
        A[N + i, (4*(i-1))+2] = xx[i+1]
        A[N + i, (4*(i-1))+3] = xx[i+1]^2
        A[N + i, (4*(i-1))+4] = xx[i+1]^3
        r[N + i] = yy[i+1]
    end
    
        # First derivative constraints: C_i'(x_i) = C_{i+1}'(x_i) or in julia notation C_i'(x_{i+1}) = C_{i+1}'(x_{i+1})
        # (note that this only goes up to N-1)
    for i in 1:N-1
        A[2*N + i, (4*(i-1))+1] = 0
        A[2*N + i, (4*(i-1))+2] = 1
        A[2*N + i, (4*(i-1))+3] = 2*xx[i+1]
        A[2*N + i, (4*(i-1))+4] = 3*xx[i+1]^2
        A[2*N + i, (4*(i-1))+5] =-A[2*N + i, (4*(i-1))+1]
        A[2*N + i, (4*(i-1))+6] =-A[2*N + i, (4*(i-1))+2]
        A[2*N + i, (4*(i-1))+7] =-A[2*N + i, (4*(i-1))+3]
        A[2*N + i, (4*(i-1))+8] =-A[2*N + i, (4*(i-1))+4]
        r[2*N + i] = 0
    end
        
        # Second derivative constraints: C_i''(x_i) = C_{i+1}''(x_i) or in julia notation C_i"(x_{i+1}) = C_{i+1}"(x_{i+1})
    for i in 1:N-1
        A[3*N + i, (4*(i-1))+1] = 0
        A[3*N + i, (4*(i-1))+2] = 0
        A[3*N + i, (4*(i-1))+3] = 2
        A[3*N + i, (4*(i-1))+4] = 6*xx[i+1]
        A[3*N + i, (4*(i-1))+5] =-A[3*N + i, (4*(i-1))+1]
        A[3*N + i, (4*(i-1))+6] =-A[3*N + i, (4*(i-1))+2]
        A[3*N + i, (4*(i-1))+7] =-A[3*N + i, (4*(i-1))+3]
        A[3*N + i, (4*(i-1))+8] =-A[3*N + i, (4*(i-1))+4]
        r[3*N + i] = 0
    end
        
        # Natural constraints for the leftmost and rightmost data point...
        # You can put these in rows 3*N and 4*N, because those aren't filled in yet!
     A[3*N, 3] = 2;
     A[3*N, 4] = 6*xx[1];
     
     A[4*N, (4*(N-1))+3] = 2;
     A[4*N, (4*(N-1))+4] = 6*xx[N+1];
    c = A \ r
    return c
end;

# Now, for Complicated Reasons, Julia can sometimes make bad decisions about how many CPU cores or threads to use.
# Turn that off for simplicity in measuring our timing.
using LinearAlgebra
LinearAlgebra.BLAS.set_num_threads(1)
     
t = @elapsed get_spline_coefficients(xx, yy)
t

# This range of matrix sizes worked okay for me (the largest one took ~5 seconds),
# but feel free to change it if it's taking too long!
sizes = range(100, 2000, 20)
times = Vector{Float64}()
for s in sizes
    s = Int(floor(s))
    xx = range(min_x, max_x, s)
    yy = sin.(xx)
    t = @elapsed get_spline_coefficients(xx, yy)
    println("Size 
t")
    push!(times, t)
end

f = Figure()
ax = Axis(f[1,1], title="Cost of spline fitting", xlabel="Size", ylabel="Compute time (sec)")
# You can make the plot a log-log plot by adding this to the "Axis" call:   yscale=log10, xscale=log10
scatter!(sizes, times)
f

# As the problem gets larger, the timing cost grows significantly faster than linearly.
# This may be due to the fact that everytime we assume bigger dimensions, we have to carry more
# zero elements in the A matrix, leading to the last line in the function (c=A/r) 
# being the most costly step. We can find the power by which the computing time grows using linear 
# regression. Specifically:

using LsqFit

# Define the power-law model
power_law(x, p) = p[1] .* sizes .^ p[2]

# Initial guess for parameters [a, b]
initial_params = [1e-5, 1e-5]

# Perform nonlinear fitting
fit = curve_fit(power_law, sizes, times, initial_params)

# Extract parameters
a, b = fit.param

println("a = $a, b = $b")

# Generate the fitted curve
fitted_y = power_law(x, fit.param)

# Plot the data and the fitted curve
Plots.scatter(sizes, times,  label="Data", xlabel="x", ylabel="y", title="Nonlinear Power-Law Fit")
Plots.plot!(sizes, fitted_y, label="Fitted Curve", linewidth=2)


# We find that: a = 3.712005966260924e-9, b = 2.8299537558847407, meaning that the power law is:
# y=a*x^b

####################################################################

using SparseArrays

# You can create a sparse matrix by passing in the ROWS, COLUMNS, and VALUES:
rows = [1, 1, 2, 3]
cols = [1, 3, 2, 3]
values = [7., 100., 23., 49.]
A = sparse(rows, cols, values)

function get_spline_coefficients_sparse(xx, yy)
    N = length(xx) - 1
    r = zeros(4*N)
    
    A_rows = Vector{Int}()
    A_cols = Vector{Int}()
    A_vals = Vector{Float64}()

    ### FILL IN THE A MATRIX HERE!
    #
    # If you used to have
    #     A[i, j] = k
    #
    # You'll now write
       
    #
    # so every line of your code before becomes three lines now -- this function will be long!
    #
    for i in 1:N  # N here corresponds to the polynomial C[i] where i goes from 1 to N 
        push!(A_rows, i)
        push!(A_cols, (4*(i-1))+1)
        push!(A_vals, 1)

        push!(A_rows, i)
        push!(A_cols, (4*(i-1))+2)
        push!(A_vals, xx[i])
       
        push!(A_rows, i)
        push!(A_cols, (4*(i-1))+3)
        push!(A_vals, xx[i]^2)
      
        push!(A_rows, i)
        push!(A_cols, (4*(i-1))+4)
        push!(A_vals, xx[i]^3)

        r[i] = yy[i]
    end
        
        # Right endpoint constraints: C_i(x_i) = y_i or in julia notation C_i(x_{i+1}) = y_{i+1}
        # Put these in the next N rows, N+1 to 2*N
    for i in 1:N
        push!(A_rows, N + i)
        push!(A_cols, (4*(i-1))+1)
        push!(A_vals, 1)
        
        push!(A_rows, N + i)
        push!(A_cols, (4*(i-1))+2)
        push!(A_vals, xx[i+1])
       
        push!(A_rows, N + i)
        push!(A_cols, (4*(i-1))+3)
        push!(A_vals, xx[i+1]^2)
      
        push!(A_rows, N + i)
        push!(A_cols, (4*(i-1))+4)
        push!(A_vals, xx[i+1]^3)

        r[N + i] = yy[i+1]
    end
    
        # First derivative constraints: C_i'(x_i) = C_{i+1}'(x_i) or in julia notation C_i'(x_{i+1}) = C_{i+1}'(x_{i+1})
        # (note that this only goes up to N-1)
    for i in 1:N-1
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+1)
        push!(A_vals, 0)
    
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+2)
        push!(A_vals, 1)
      
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+3)
        push!(A_vals, 2*xx[i+1])
       
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+4)
        push!(A_vals, 3*xx[i+1]^2)
    
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+5)
        push!(A_vals, 0)
        
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+6)
        push!(A_vals,-1)
       
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+7)
        push!(A_vals,-2*xx[i+1])
        
        push!(A_rows, 2*N + i)
        push!(A_cols, (4*(i-1))+8)
        push!(A_vals,-3*xx[i+1]^2)

        r[2*N + i] = 0
    end
        
        # Second derivative constraints: C_i''(x_i) = C_{i+1}''(x_i) or in julia notation C_i"(x_{i+1}) = C_{i+1}"(x_{i+1})
    for i in 1:N-1
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+1)
        push!(A_vals, 0)
       
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+2)
        push!(A_vals, 0)
       
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+3)
        push!(A_vals, 2)
       
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+4)
        push!(A_vals, 6*xx[i+1])
      
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+5)
        push!(A_vals, 0)

        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+6)
        push!(A_vals, 0)
    
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+7)
        push!(A_vals,-2)
       
        push!(A_rows, 3*N + i)
        push!(A_cols, (4*(i-1))+8)
        push!(A_vals,-6*xx[i+1])

        r[3*N + i] = 0
    end
        
        # Natural constraints for the leftmost and rightmost data point...
        # You can put these in rows 3*N and 4*N, because those aren't filled in yet!
     push!(A_rows, 3*N)
     push!(A_cols, 3)
     push!(A_vals, 2)

     push!(A_rows, 3*N)
     push!(A_cols, 4)
     push!(A_vals, 6*xx[1])
     
     push!(A_rows, 4*N)
     push!(A_cols, (4*(N-1))+3)
     push!(A_vals, 2)
   
     push!(A_rows, 4*N)
     push!(A_cols, (4*(N-1))+4)
     push!(A_vals, 6*xx[N+1])
    
    A = sparse(A_rows, A_cols, A_vals)
    c = A \ r
    return c
end

c_sparse = get_spline_coefficients_sparse(xx, yy)
spline_sparse = make_spline(xx, c_sparse);

# Now go ahead and plot your resulting spline!
f = Figure()
ax = Axis(f[1,1], title="Data and Spline", xlabel="X", ylabel="Y")
scatter!(xx, yy)
plotx = range(min_x, max_x, 1000)
lines!(plotx, spline, linewidth=3)
lines!(plotx, spline_sparse)
f

sparse_sizes = range(100, 2000, 20)
sparse_times = Vector{Float64}()
for s in sparse_sizes
    s = Int(floor(s))
    xx = range(min_x, max_x, s)
    yy = sin.(xx)
    t = @elapsed get_spline_coefficients_sparse(xx, yy)
    println("Size 
t")
    push!(sparse_times, t)
end

f = Figure()
ax = Axis(f[1,1], title="Cost of spline fitting", xlabel="Size", ylabel="Compute time (sec)",
        yscale=log10)
scatter!(sizes, times)
scatter!(sparse_sizes, sparse_times)
f

# The sparse version is significantly faster. The curve looks like a power law again
# but with smaller power a. We perform the same power law fitting as before and get the 
# results shown below

using LsqFit

# Define the power-law model
power_law(x, p) = p[1] .* sparse_sizes .^ p[2]

# Perform nonlinear fitting
fit_ = curve_fit(power_law, sparse_sizes, sparse_times, initial_params)

# Extract parameters
a_, b_ =  fit_.param

println("a_ = $a_, b_ = $b_")

# Generate the fitted curve
fitted_y_= power_law(x, fit_.param)

# Plot the data and the fitted curve

Plots.scatter(sizes, log.(times),  label="Data", xlabel="x", ylabel="y", title="Nonlinear Power-Law Fit")
Plots.plot!(sizes, log.(fitted_y), label="Fitted Curve", linewidth=2)
Plots.scatter!(sparse_sizes, log.(sparse_times), label="Data Sparse")
Plots.plot!(sparse_sizes, log.(fitted_y_), label="Fitted Curve Sparse", linewidth=2)

# So the power law coefficients for th sparse version are: a_ = 2.0348783789076636e-5, b_ = 0.8088423989978681
# meaning that the spaese version did indeed reduce the cost of computing the spline but also reduced the rate 
# at which this cost increases by using bigger dimensions.  
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

f = Figure()
Axis(f[1,1], title="Data and Spline", xlabel="X", ylabel="Y")
scatter!(xx, yy)
plotx = range(min_x, max_x, 1000)
lines!(plotx, spl(plotx), linewidth=3)
lines!(plotx, sin.(plotx))
Axis(f[2,1], title="Error in spline", xlabel="X", ylabel="delta-Y")
lines!(plotx, spl(plotx) .- sin.(plotx))
scatter!(xx, zeros(length(xx)))
f

N = length(xx) - 1
A = zeros((4*N, 4*N))
r = zeros(4*N);

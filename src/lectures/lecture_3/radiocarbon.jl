# ## A simple model for $R \equiv ^{14}\!\!C/^{12}C$ normalized so that $R=1$ in the atmosphere
#
# $$ \frac{\partial R}{\partial t} + [\mathbf{TRdiv}]R = \kappa\boldsymbol{\Lambda}(1-R) -\lambda R $$

# Here we will redo the radiocarbon age example, but using the OCIM transport operator made available through the AIBECS package
using AIBECS
mask, grd, T_OCIM = OCIM1.load() ;

# Some useful OCIM stuff
iwet = findall(x -> x == 1, vec(mask));  # index to wet gridboxes
dz1 = grd["dzt"][1]                      # thickness of the top layer
z = vec(grd["ZT3d"])[iwet]               # depth of the gridbox centers
nwet = length(iwet)                      # number of wet grid boxes
dv = vec(grd["DZT3d"])[iwet].*vec(grd["DYT3d"])[iwet].*vec(grd["DXT3d"])[iwet]; # volume of the gridboxes
lat, lon = vec(grd["yt"]), vec(grd["xt"]); # latitudes and longitudes (useful for plotting)
depth = vec(grd["zt"]) ;


# Make a table of parameters
t = empty_parameter_table()               # initialize table of parameters
add_parameter!(t, :λ, 1 / (5730*log(2))u"yr") # add the radioactive decay e-folding timescale
add_parameter!(t, :κ, 50u"m" / 10u"yr")
initialize_Parameters_type(t, "C14_Parameters")             # Generate the parameter table

# The model parameters
p₀ = C14_Parameters()

# The sources and sinks of $R$
function sms_14c(R, p) 
    λ = p.λ
    κ = p.κ
    return κ * (z .< 20) .* (1.0 .- R) / dz1 - λ.*R
end
# The tracer transport operator
T_14c(p) = T_OCIM


# Prepare the stuff needed for the AIBECS solver:
T_matrices = (T_14c,)            # bundles all the transport matrices in a tuple
sources_minus_sinks = (sms_14c,) # bundles all the source-sink functions in a tuple
F, ∇ₓF = state_function_and_Jacobian(T_matrices, sources_minus_sinks, nwet) # generates the state function (and its Jacobian!)
#
x₀ = ones(nwet)                           # initial iterate for the solver
prob = SteadyStateProblem(F, ∇ₓF, x₀, p₀) # define the problem
R = solve(prob, CTKAlg())                 # solve the problem
c14age = -log.(R)/p₀.λ;                   # convert R to $^{14}C-age$


# Make some plots
c14age_3d = NaN * mask     # creates a 3D array of NaNs
c14age_3d[iwet] = c14age   # Fills the wet grid boxes with the age values
size(c14age_3d)            # Just to check the size of age_3D

# Pick a layer to plot 
iz = findfirst(depth .> 700) # aim for a depth of ~ 700 m
iz, depth[iz]
c14age_3d_1000m_yr = c14age_3d[:,:,iz] * ustrip(1.0u"s" |> u"yr")
#

#
ENV["MPLBACKEND"]="qt5agg"
using PyPlot, PyCall
#
clf()
ccrs = pyimport("cartopy.crs")
ax = subplot(projection=ccrs.Robinson(central_longitude=-155.0))
ax.coastlines()
lon_cyc = [lon; 360+lon[1]] # making it cyclic for Cartopy
age_cyc = hcat(c14age_3d_1000m_yr, c14age_3d_1000m_yr[:,1])
p = contourf(lon_cyc, lat, age_cyc, levels=0:100:1800, transform=ccrs.PlateCarree(), zorder=-1)
colorbar(p, orientation="horizontal")
gcf() # gets the current figure to display

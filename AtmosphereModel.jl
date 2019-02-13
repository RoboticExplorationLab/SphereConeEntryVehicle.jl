using ApproxFun
using LinearAlgebra
using PyCall
using PyPlot
using JLD

function compute_chebyshev_coefficients_atmosphere_log(height, density, order)
    altitude_min = height[1]
    altitude_max = height[end]
    num_nodes = size(density)[1]
    space = Chebyshev(Interval(altitude_min, altitude_max))
    # a non-default grid
    points = range(altitude_min, stop=altitude_max, length=num_nodes) 
    values = log.(density)
    # Create a Vandermonde matrix by evaluating the basis at the grid
    V = zeros(Float64, num_nodes, order)
    for k = 1:order
        V[:, k] = Fun(space,[zeros(k-1);1]).(points)
    end
    atmosphere_polynomial_log = Fun(space, V\values)
    atmosphere_coefficients = atmosphere_polynomial_log.coefficients
    return atmosphere_coefficients
end

function save_chebyshev_coefficients_atmosphere_log(height, density, order, filename)
    atmosphere_coefficients = compute_chebyshev_coefficients_atmosphere_log(height, density, order)
    save(filename, "atmosphere_coefficients", atmosphere_coefficients, 
    	"altitude_min", altitude_min, "altitude_max", altitude_max)
end

function load_chebyshev_approximation_atmosphere_log(filename)
    file = load(filename)
#     if altitude_min != file["altitude_min"] || altitude_max != file["altitude_max"]
#         println("The altitude parameters (min/max) do not correspond to the loaded file.")
#     end
    altitude_min = file["altitude_min"]
    altitude_max = file["altitude_max"]
    atmosphere_coefficients = file["atmosphere_coefficients"]
    return atmosphere_coefficients, altitude_min, altitude_max
end

function atmosphere_density_chebyshev(altitude)
    atmosphere_coefficients, altitude_min, altitude_max = 
    	load_chebyshev_approximation_atmosphere_log("chebyshev_coefficients_atmosphere.jld")
    space = Chebyshev(Interval(altitude_min, altitude_max))
    density_polynomial_log = Fun(space, atmosphere_coefficients)
    density = exp(density_polynomial_log(altitude))
    return density
end

function compute_atmosphere_density_exp(altitude)
    # Exponential model
    ρ_inf = ρ_0 * exp(- altitude / h)
    return ρ_inf
end

function compute_atmosphere_density_parametric(altitude)
    temperature = -23.4 - 0.00222 * altitude
    pressure = 0.699 * exp(-0.00009 * altitude)
    ρ_inf = pressure / (0.1921 * (temperature + 273.1))
    return ρ_inf
end

function plot_atmosphere_density_empirical(density, height, num_nodes; display=true)
    altitude_min = height[1]
    altitude_max = height[end]    
    height_poly = range(altitude_min, stop=altitude_max, length=num_nodes)
    density_poly = atmosphere_density_chebyshev.(height_poly)
    plot(height_poly, log.(10, density_poly), color="blue", linewidth=3.0, 
    	linestyle="-.", label="density_poly")
    plot(height, log.(10, density), color="red", linewidth=2.0, 
    	linestyle="-", label="density_table")
    legend()
    xlabel("Altitude in m")
    ylabel("log_10(density) in kg/m^3")
    title("Mars Atmosphere Density.")
    grid("on")
    if display 
        tight_layout() 
        show()
    end
    return
end

function plot_atmosphere_density_comparison(altitude_min, altitude_max, num_nodes; display=true)
    altitude = range(altitude_min, stop=altitude_max, length=num_nodes)

    density_exp = compute_atmosphere_density_exp.(altitude)
    density_parametric = compute_atmosphere_density_parametric.(altitude)
    density_chebyshev = atmosphere_density_chebyshev.(altitude)

    plot(altitude, log.(10, density_exp), color="blue", linewidth=2.0, 
    	linestyle="-", label="density_exp")
    plot(altitude, log.(10, density_parametric), color="red", linewidth=2.0, 
    	linestyle=":", label="density_parametric")
    plot(altitude, log.(10, density_chebyshev), color="green", linewidth=2.0, 
    	linestyle="-.", label="density_chebyshev")
    legend()
    xlabel("Altitude in m")
    ylabel("log_10(density) in kg/m^3")
    title("Mars Atmosphere Density.")
    grid("on")
    if display
        tight_layout() 
        show()
    end
    return
end

if PROGRAM_FILE == @__FILE__ 
	# Data used for the three models
	# Extracted from https://ntrs.nasa.gov/archive/nasa/casi.ntrs.nasa.gov/20070032693.pdf.
	altitude_min = 0
	altitude_max = 150000
	num_nodes = 16
	height = range(altitude_min, stop=altitude_max, length=num_nodes)
	density = [1.550e-2
	            6.470e-3
	            2.630e-3
	            9.800e-4
	            3.400e-4
	            1.080e-4
	            3.180e-5
	            8.730e-6
	            2.290e-6
	            6.010e-7
	            1.590e-7
	            4.140e-8
	            1.190e-8
	            3.760e-9
	            1.090e-9
	            4.730e-10]
	h = 10800 # Atmospheric Scale Height [m]
	ρ_0 = 0.020 # Surface density of Mars atmosphere [kg.m^-3]

	# Saving the polynomial model. 
	order = 10
	filename = "chebyshev_coefficients_atmosphere.jld"
	save_chebyshev_coefficients_atmosphere_log(height, density, order, filename)

	#Plotting the curves
	altitude_min = 5000
	altitude_max = 112000
	num_nodes_plot = 1000
	plot_atmosphere_density_empirical(density, height, num_nodes)
	plot_atmosphere_density_comparison(altitude_min, altitude_max, num_nodes_plot)
end

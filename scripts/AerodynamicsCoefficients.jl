using ApproxFun
using LinearAlgebra
using HCubature
using StaticArrays
using JLD

function C_X2_sphere_parameterized_integrand(x::SVector{2,Float64}, r_d, α)
	u = x[1]
	v = x[2]
	z = (u * cos(α) + sin(α) * sqrt(r_d^2 - u^2) * cos(v)) / r_d
	postive_part_z = max(0, z)
	y = postive_part_z^2 * u
	return y
end

function C_Y2_sphere_parameterized_integrand(x::SVector{2,Float64}, r_d, α)
	u = x[1]
	v = x[2]
	z = (u * cos(α) + sin(α) * sqrt(r_d^2 - u^2) * cos(v)) / r_d
	postive_part_z = max(0, z)
	y = postive_part_z^2 * cos(v) * sqrt(r_d^2 - u^2)
	return y
end

function C_X2_cone_parameterized_integrand(x::SVector{2,Float64}, δ, α)
	u = x[1]
	v = x[2]
	z = (cos(α) + sin(α) / tan(δ) * cos(v)) / sqrt(1 + 1.0 / tan(δ)^2)
	postive_part_z = max(0, z)
	y = postive_part_z^2 * u
	return y
end

function C_Y2_cone_parameterized_integrand(x::SVector{2,Float64}, δ, α)
	u = x[1]
	v = x[2]
	z = (cos(α) + sin(α) / tan(δ) * cos(v)) / sqrt(1 + 1.0 / tan(δ)^2)
	postive_part_z = max(0, z)
	y = postive_part_z^2 * u * cos(v) / tan(δ)
	return y
end

function C_n2_cone_parameterized_integrand(x::SVector{2,Float64}, δ, α, l)
	u = x[1]
	v = x[2]
	z = (cos(α) + sin(α) / tan(δ) * cos(v)) / sqrt(1 + 1.0 / tan(δ)^2)
	postive_part_z = max(0, z)
	second_term = u * cos(v) * ((l - u / tan(δ)) / tan(δ) - u)
	y = postive_part_z^2 * second_term
	return y
end

function aerodynamics_coefficents(α, r_max, r_min, l, x_g)
	# This function return the aerodynamics coefficents associated to forces and moments for the 
	# sphere-cone geometry. Important to notice that the dimension of the output are not as defined
	# in the original paper. 
	# C_F2 [m^2] 
	# C_τ2 [m^3]
	# In the paper, 
	# C_F2 * A_ref [m^2] 
	# C_τ2 * A_ref * l_ref[m^3]
	# We fuse the coefficent with the area and length of reference since they are not useful for 
	
	# Geometry of the vehicle
	# angle of the cone
	δ = atan(r_max, l)
	# radius of the sphere
	r_d = r_min / cos(δ)
	# distance between the center of the sphere and the base of the small cone
	r_0 = sqrt(r_d^2 - r_min^2)
	# length of the section of cone.
	f = l * (1 - r_min/r_max)
	# vector from center of mass to center of the sphere expressed in freestream frame
	GS = [f - r_0 - x_g , 0, 0]
	# vector from center of mass to center of large cone base expressed in freestream frame
	GC = [- x_g, 0, 0]

	C_F2_cone, C_τ2_cone = aerodynamics_cone_coefficents(α, r_max, r_min, l)
	C_F2_sphere, C_τ2_sphere = aerodynamics_sphere_coefficents(α, r_max, r_min, l)

	C_F2 = C_F2_cone + C_F2_sphere
	C_τ2 = [0, 0, C_τ2_cone[3] + C_τ2_sphere[3]] + cross(GC, C_F2_cone) + cross(GS, C_F2_sphere)
	return C_F2, C_τ2
end

function aerodynamics_cone_coefficents(α, r_max, r_min, l)
	# This function return the aerodynamics coefficents associated to forces and moments for the 
	# cone geometry. Important to notice that the dimension of the output are not as defined
	# in the original paper. 
	# C_F2 [m^2] 
	# C_τ2 [m^3]
	# In the paper, 
	# C_F2 * A_ref [m^2] 
	# C_τ2 * A_ref * l_ref[m^3]
	# We fuse the coefficent with the area and length of reference since they are not useful for 
	# forces and moments computation. 

	# Geometry of the vehicle
	# angle of the cone
	δ = atan(r_max, l)
	# radius of the sphere
	r_d = r_min / cos(δ)
	# distance between the center of the sphere and the base of the small cone
	r_0 = sqrt(r_d^2 - r_min^2)
	# length of the section of cone.
	f = l * (1 - r_min/r_max)

	# Aerodynamics coefficients for the cone 
	start_bound = [r_min, 0.0]
	end_bound = [r_max, 2 * pi]

	function C_X2_cone_integrand(x::SVector{2,Float64})
		return C_X2_cone_parameterized_integrand(x, δ, α)
	end
	function C_Y2_cone_integrand(x::SVector{2,Float64})
		return C_Y2_cone_parameterized_integrand(x, δ, α)
	end
	function C_n2_cone_integrand(x::SVector{2,Float64})
		return C_n2_cone_parameterized_integrand(x, δ, α, l)
	end

	# Parameters of integration to limit memory usage and running time.
	rtol = 1e-4
	maxevals = 1000

	# Integration
	integral_X, error_X = hcubature(C_X2_cone_integrand, start_bound, end_bound;
		norm=norm, rtol=rtol, maxevals=maxevals)
	C_X2 = - 2.0 * integral_X
	integral_Y, error_Y = hcubature(C_Y2_cone_integrand, start_bound, end_bound;
		norm=norm, rtol=rtol, maxevals=maxevals)
	C_Y2 = - 2.0 * integral_Y
	integral_n, error_n = hcubature(C_n2_cone_integrand, start_bound, end_bound;
		norm=norm, rtol=rtol, maxevals=maxevals)
	C_n2 = - 2.0 * integral_n

	C_F2 = [C_X2, C_Y2, 0]
	C_τ2 = [0, 0, C_n2]
	return C_F2, C_τ2
end

function aerodynamics_sphere_coefficents(α, r_max, r_min, l)
	# This function return the aerodynamics coefficents associated to forces and moments for the 
	# sphere geometry. Important to notice that the dimension of the output are not as defined
	# in the original paper. 
	# C_F2 [m^2] 
	# C_τ2 [m^3]
	# In the paper, 
	# C_F2 * A_ref [m^2] 
	# C_τ2 * A_ref * l_ref[m^3]
	# We fuse the coefficent with the area  and length of reference since they are not useful for 
	# forces and moments computation. 

	# Geometry of the vehicle
	# angle of the cone
	δ = atan(r_max, l)
	# radius of the sphere
	r_d = r_min / cos(δ)
	# distance between the center of the sphere and the base of the small cone
	r_0 = sqrt(r_d^2 - r_min^2)
	# length of the section of cone.
	f = l * (1 - r_min/r_max)

	# define domain of integration and number of evaluation in the integral computation.
	start_bound = [r_0, 0.0]
	end_bound = [r_d, 2 * pi]
	# Define both integrands with fixed parameters r_d and α.
	function C_X2_sphere_integrand(x::SVector{2,Float64})
		return C_X2_sphere_parameterized_integrand(x, r_d, α)
	end
	function C_Y2_sphere_integrand(x::SVector{2,Float64})
		return C_Y2_sphere_parameterized_integrand(x, r_d, α)
	end

	# Parameters of integration to limit memory usage and running time.
	rtol = 1e-4
	maxevals = 1000
	# Integration
	integral_X, error_X = hcubature(C_X2_sphere_integrand, start_bound, end_bound;
		norm=norm, rtol=rtol, maxevals=maxevals)
	C_X2 = - 2.0 * integral_X
	integral_Y, error_Y = hcubature(C_Y2_sphere_integrand, start_bound, end_bound;
		norm=norm, rtol=rtol, maxevals=maxevals)
	C_Y2 = - 2.0 * integral_Y
	C_n2 = 0.0

	C_F2 = [C_X2, C_Y2, 0]
	C_τ2 = [0, 0, C_n2]
	return C_F2, C_τ2
end

function plot_aerodynamics_ceofficients(r_max, r_min, l, x_g, num_nodes; display=true)
	α_vect = range(0, stop=pi/2, length=num_nodes)
	C_F2_vect = zeros(num_nodes, 3)
	C_τ2_vect = zeros(num_nodes, 3)
	for i=1:num_nodes
		C_F2_vect[i, :], C_τ2_vect[i, :] = aerodynamics_coefficents(α_vect[i], r_max, r_min, l, x_g)
	end
	plot(α_vect, C_F2_vect[:, 1], color="blue", linewidth=2.0, linestyle="-", label="C_X2_vehicle")
	plot(α_vect, C_F2_vect[:, 2], color="red", linewidth=2.0, linestyle="-", label="C_Y2_vehicle")
	plot(α_vect, C_τ2_vect[:, 3], color="green", linewidth=2.0, linestyle="-", label="C_n2_vehicle")
    plot(α_vect, C_F2_vect[:, 2] ./ tan.(α_vect), 
    	color="red", linewidth=2.0, linestyle=":", label="C_Y2_vehicle / tan(α)")
    plot(α_vect, C_τ2_vect[:, 3] ./ tan.(α_vect), 
    	color="green", linewidth=2.0, linestyle=":", label="C_n2_vehicle / tan(α)")
	legend()
	xlabel("α")
	title("Aerodynamics coefficients for the vehicle.")
	grid("on")
	if display 
		tight_layout()
		show()
	end
	return
end

function plot_aerodynamics_ceofficients_cone(r_max, r_min, l, num_nodes; display=true)
	α_vect = range(0, stop=pi/2, length=num_nodes)
	C_F2_vect = zeros(num_nodes, 3)
	C_τ2_vect = zeros(num_nodes, 3)
	for i=1:num_nodes
		C_F2_vect[i, :], C_τ2_vect[i, :] = aerodynamics_cone_coefficents(α_vect[i], r_max, r_min, l)
	end
	plot(α_vect, C_F2_vect[:, 1], color="blue", linewidth=2.0, linestyle="-", label="C_X2_cone")
	plot(α_vect, C_F2_vect[:, 2], color="red", linewidth=2.0, linestyle="-", label="C_Y2_cone")
	plot(α_vect, C_τ2_vect[:, 3], color="green", linewidth=2.0, linestyle="-", label="C_n2_cone")
    plot(α_vect, C_F2_vect[:, 2] ./ tan.(α_vect), 
    	color="red", linewidth=2.0, linestyle=":", label="C_Y2_cone / tan(α)")
    plot(α_vect, C_τ2_vect[:, 3] ./ tan.(α_vect), 
    	color="green", linewidth=2.0, linestyle=":", label="C_n2_cone / tan(α)")
	legend()
	xlabel("α")
	title("Aerodynamics coefficients for the cone.")
	grid("on")
	if display
		tight_layout() 
		show()
	end
	return
end

function plot_aerodynamics_ceofficients_sphere(r_max, r_min, l, num_nodes; display=true)
	α_vect = range(0, stop=pi/2, length=num_nodes)
	C_F2_vect = zeros(num_nodes, 3)
	C_τ2_vect = zeros(num_nodes, 3)
	for i=1:num_nodes
		C_F2_vect[i, :], C_τ2_vect[i, :] = 
			aerodynamics_sphere_coefficents(α_vect[i], r_max, r_min, l)
	end
	plot(α_vect, C_F2_vect[:, 1], color="blue", linewidth=2.0, linestyle="-", label="C_X2_sphere")
	plot(α_vect, C_F2_vect[:, 2], color="red", linewidth=2.0, linestyle="-", label="C_Y2_sphere")
	plot(α_vect, C_τ2_vect[:, 3], color="green", linewidth=2.0, linestyle="-", label="C_n2_sphere")
	plot(α_vect, C_F2_vect[:, 2] ./ tan.(α_vect), 
		color="red", linewidth=2.0, linestyle=":", label="C_Y2_sphere / tan(α)")
	plot(α_vect, C_τ2_vect[:, 3] ./ tan.(α_vect), 
		color="green", linewidth=2.0, linestyle=":", label="C_n2_sphere / tan(α)")
	legend()
	xlabel("α")
	title("Aerodynamics coefficients for the sphere.")
	grid("on")
	if display 
		tight_layout()
		show()
	end
	return
end

function compute_chebyshev_coefficients(r_max, r_min, l, x_g, order)
	space = Chebyshev(Interval(-pi/2, pi/2))
	p = points(space, order) # the default grid
	function compute_C_X2(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_X2 = C_F2[1]
		return C_X2
	end
	function compute_C_Y2(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_Y2 = C_F2[2]
		return C_Y2
	end
	function compute_C_n2(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_n2 = C_τ2[3]
		return C_n2
	end

	C_X2_values = compute_C_X2.(p)      # values at the default grid
	C_Y2_values = compute_C_Y2.(p)      # values at the default grid
	C_n2_values = compute_C_n2.(p)      # values at the default grid
	C_X2_polynomial = Fun(space, ApproxFun.transform(space, C_X2_values))
	C_Y2_polynomial = Fun(space, ApproxFun.transform(space, C_Y2_values))
	C_n2_polynomial = Fun(space, ApproxFun.transform(space, C_n2_values))

	coefficients = zeros(order, 3)
	coefficients[:, 1] =  C_X2_polynomial.coefficients
	coefficients[:, 2] =  C_Y2_polynomial.coefficients
	coefficients[:, 3] =  C_n2_polynomial.coefficients
	return coefficients
end

function compute_chebyshev_coefficients_tan(r_max, r_min, l, x_g, order)
	space = Chebyshev(Interval(-pi/2, pi/2))
	p = points(space, order) # the default grid
	function compute_C_X2(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_X2 = C_F2[1]
		return C_X2
	end
	function compute_C_Y2_tan(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_Y2 = C_F2[2]
		return C_Y2 / tan(α)
	end
	function compute_C_n2_tan(α)
		C_F2, C_τ2 = aerodynamics_coefficents(α, r_max, r_min, l, x_g)
		C_n2 = C_τ2[3]
		return C_n2 / tan(α)
	end

	C_X2_values = compute_C_X2.(p)      # values at the default grid
	C_Y2_tan_values = compute_C_Y2_tan.(p)      # values at the default grid
	C_n2_tan_values = compute_C_n2_tan.(p)      # values at the default grid
	C_X2_polynomial = Fun(space, ApproxFun.transform(space, C_X2_values))
	C_Y2_tan_polynomial = Fun(space, ApproxFun.transform(space, C_Y2_tan_values))
	C_n2_tan_polynomial = Fun(space, ApproxFun.transform(space, C_n2_tan_values))

	coefficients_tan = zeros(order, 3)
	coefficients_tan[:, 1] =  C_X2_polynomial.coefficients
	coefficients_tan[:, 2] =  C_Y2_tan_polynomial.coefficients
	coefficients_tan[:, 3] =  C_n2_tan_polynomial.coefficients
	return coefficients_tan
end

function save_chebyshev_coefficients(r_max, r_min, l, x_g, order, filename)
	coefficients = compute_chebyshev_coefficients(r_max, r_min, l, x_g, order)
	save(filename, "coefficients", coefficients)
end

function save_chebyshev_coefficients_tan(r_max, r_min, l, x_g, order, filename)
	coefficients_tan = compute_chebyshev_coefficients_tan(r_max, r_min, l, x_g, order)
	save(filename, "r_max", r_max, "r_min", r_min, "l", l, "x_g", x_g, "order", order,
		"coefficients_tan", coefficients_tan)
end

function load_chebyshev_approximation(filename)
	file = load(filename)
	if r_max != file["r_max"] || r_min != file["r_min"] || l != file["l"] || x_g != file["x_g"]
		println("The geometric parameters does not correspond to the loaded file.")
		return
	end
	coefficients = file["coefficients"]
	return coefficients
end

function load_chebyshev_approximation_tan(filename, r_max, r_min, l, x_g)
	file = load(filename)
	if r_max != file["r_max"] || r_min != file["r_min"] || l != file["l"] || x_g != file["x_g"]
		println("The geometric parameters do not correspond to the loaded file.")
		return
	end
	coefficients_tan = file["coefficients_tan"]
	return coefficients_tan
end

function aerodynamics_coefficents_chebyshev(α)
	coefficients = load_chebyshev_approximation("chebyshev_coefficients.jld")
	C_X2_coefficients = coefficients[:, 1]
	C_Y2_coefficients = coefficients[:, 2]
	C_n2_coefficients = coefficients[:, 3]

	order = size(coefficients)[1]
	space = Chebyshev(Interval(-pi/2, pi/2))
	C_X2 = Fun(space, C_X2_coefficients)
	C_Y2 = Fun(space, C_Y2_coefficients)
	C_n2 = Fun(space, C_n2_coefficients)

	C_F2 = [C_X2(α), C_Y2(α), 0]
	C_τ2 = [0, 0, C_n2(α)]
	return C_F2, C_τ2
end

function aerodynamics_coefficents_chebyshev_tan(α, r_max, r_min, l, x_g)
	coefficients_tan = load_chebyshev_approximation_tan(
		"chebyshev_coefficients.jld", r_max, r_min, l, x_g)
	C_X2_coefficients = coefficients_tan[:, 1]
	C_Y2_tan_coefficients = coefficients_tan[:, 2]
	C_n2_tan_coefficients = coefficients_tan[:, 3]

	order = size(coefficients_tan)[1]
	space = Chebyshev(Interval(-pi/2, pi/2))
	C_X2 = Fun(space, C_X2_coefficients)
	C_Y2_tan = Fun(space, C_Y2_tan_coefficients)
	C_n2_tan = Fun(space, C_n2_tan_coefficients)

	C_F2_tan = [C_X2(α), C_Y2_tan(α), 0]
	C_τ2_tan = [0, 0, C_n2_tan(α)]
	return C_F2_tan, C_τ2_tan
end

function aerodynamics_functions_chebyshev_tan(r_max, r_min, l, x_g)
	coefficients_tan = load_chebyshev_approximation_tan(
		"chebyshev_coefficients.jld", r_max, r_min, l, x_g)
	C_X2_coefficients = coefficients_tan[:, 1]
	C_Y2_tan_coefficients = coefficients_tan[:, 2]
	C_n2_tan_coefficients = coefficients_tan[:, 3]

	order = size(coefficients_tan)[1]
	space = Chebyshev(Interval(-pi/2, pi/2))
	C_X2 = Fun(space, C_X2_coefficients)
	C_Y2_tan = Fun(space, C_Y2_tan_coefficients)
	C_n2_tan = Fun(space, C_n2_tan_coefficients)
	return C_X2, C_Y2_tan, C_n2_tan
end

# if PROGRAM_FILE == @__FILE__ 
# 	# Parameters
# 	r_min = 0.20 # smallest radius of the cone [m]
# 	r_max = 1.30 # largest radius of the cone [m]
# 	l = 0.55 # length of the cone [m]
# 	m = 569.7 # mass of the Phoenix entry system [kg]
# 	x_g = -0.152 # axial center-of-gravity location [m]
# 	order = 10
# 	filename = "chebyshev_coefficients.jld"
# 	save_chebyshev_coefficients_tan(r_max, r_min, l, x_g, order, filename)
# 	coefficients_tan = load_chebyshev_approximation_tan(filename, r_max, r_min, l, x_g)
# 	α = 0.0
# 	aerodynamics_coefficents_chebyshev_tan(α, r_max, r_min, l, x_g);
# end
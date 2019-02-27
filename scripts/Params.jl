function add_l_to_params!(params)
	params["l"] = params["r_max"] / tan(params["δ"])
	return params
end

function add_aerodynamics_to_params!(params)
	F_hat, τ_hat = illuminated_integral_aerodynamics(params["r_max"],
		params["r_min"], params["l"], params["x_g"])
	params["F_hat"] = F_hat
	params["τ_hat"] = τ_hat
	return params
end

function complete_params!(params)
	add_l_to_params!(params)
	add_aerodynamics_to_params!(params)
	return params
end

# params = Dict("r_min" => 0.20, # smallest radius of the cone [m]
# 			  "r_max" => 1.30, # largest radius of the cone [m]
# 			  "δ" => 70 / 360 * 2 * pi, # opening amgle of the cone [rad]
# 			  "m" => 569.7, # mass of the Phoenix entry system [kg]
# 			  "x_g" => -0.152, # axial center-of-gravity location [m]
#             "c" => 1.20, # distance between the radial axis and the cluster of thrusters[m]
#             "F_max" => 36.0, # maximum force applied by the thrusters [N]
#			  "Jxx" => 293.15, # Phoenix entry system [kg.m^2]
# 			  "Jyy" => 184, # Phoenix entry system [kg.m^2]
# 			  "Jzz" => 208.02, # Phoenix entry system [kg.m^2]
# 			  "Jxy" => 0.451, # Phoenix entry system [kg.m^2]
# 			  "Jxz" => -4.424, # Phoenix entry system [kg.m^2] 
# 			  "Jyz" => 0.372, # Phoenix entry system [kg.m^2]
# 			  "g" => 3.711, # Mars gravity [m^2.s^-1]
# 			  "h" => 10800, # Atmospheric Scale Height [m]
# 			  "ρ_0" => 0.020, # Surface density of Mars atmosphere [kg.m^-3]
# 			  "r_p" => 3389.5e3, # Volumetric mean radius of Mars [m]
# 			  "ω_p" => [0, 0, 7.088e-05]); # Angular velocity of Mars [rad.s^-1]

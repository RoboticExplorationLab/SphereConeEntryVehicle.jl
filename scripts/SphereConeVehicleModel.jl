using LinearAlgebra
using StaticArrays
include("AerodynamicsModel.jl")
include("ControlModel.jl")

## Utilities
"""
@(SIGNATURES)
	Conjuguate of a quaternion. 
"""
function qconj(q)
    p = [q[1]; -q[2:4]]
return p / norm(q)^2
end

"""
@(SIGNATURES)
	Rotate a vector by a quaternion
"""
function qrot(q, r)
    # with quaternion q_1/0 
    # qrot(q_1/0, r_1) = r_0
    # qrot(qconj(q_1/0), r_0) = r_1
    r + 2*cross(q[2:4],cross(q[2:4],r) + q[1]*r)
end

"""
@(SIGNATURES)
	Multiplication of two quaternions (q = [s;v])
"""
function qmult(q1, q2)
	[q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end

"""
@(SIGNATURES)
	Transformation from modified Rodrigues parameters (mrp) to rotation matrix. 
"""
function mrptodcm(mrp)
    # Converts a vector of Modified Rodrigues Parameters to a Rotation Matrix
    mrp2 = mrp' * mrp;
    S = [0       -mrp[3] mrp[2];
         mrp[3]  0       -mrp[1];
         -mrp[2] mrp[1]  0]
    I = diagm(0=>fill(1., 3))
    R = I + (8 * S * S + 4 * (1 - mrp2) * S) / ((1 + mrp2)^2)
    println(R'R)
    return R
end

function entry_vehicle_dynamics!(ẋ, X, u, params)
	## States: X ∈ R^13; q = [s;v]
	# x
	# y
	# z
	# q0
	# q1
	# q2
	# q3
	# xdot
	# ydot
	# zdot
	# omega1
	# omega2
	# omega3

	x = X[1:3]
	q = X[4:7]./norm(X[4:7]) #normalize quaternion
	v = X[8:10]
	ω = X[11:13]

	# Parameters
	r_min = params["r_min"] # smallest radius of the cone [m]
	r_max = params["r_max"] # largest radius of the cone [m]
	l = params["l"] # largest radius of the cone [m]
	m = params["m"] # mass of the Phoenix entry system [kg]
	x_g = params["x_g"] # axial center-of-gravity location [m]
	Jxx = params["Jxx"] # Phoenix entry system [kg.m^2]
	Jyy = params["Jyy"] # Phoenix entry system [kg.m^2]
	Jzz = params["Jzz"] # Phoenix entry system [kg.m^2]
	Jxy = params["Jxy"] # Phoenix entry system [kg.m^2]
	Jxz = params["Jxz"] # Phoenix entry system [kg.m^2]
	Jyz = params["Jyz"] # Phoenix entry system [kg.m^2]
	g = params["g"] # Mars gravity [m^2.s^-1]
	h = params["h"] # Atmospheric Scale Height [m]
	ρ_0 = params["ρ_0"] # Surface density of Mars atmosphere [kg.m^-3]
	r_p = params["r_p"] # Volumetric mean radius of Mars [m]
	ω_p = params["ω_p"] # Angular velocity of Mars [rad.s^-1]

	J = Matrix([Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz])# Inertia matrix [kg.m^2]
	Jinv = inv(J) # inverted inertia matrix [kg^-1.m^-2]

	# Mars atmosphere model
	altitude = norm(x) - r_p
	ρ_inf = ρ_0 * exp(- altitude / h)

	# Freestream Velocity
	# Static fluid assumption
	V_inf = - qrot(qconj(q), v + cross(x, ω_p))

	# Compute the force and torque applied of the vehicle at the center of mass.
	F_hat = params["F_hat"]
	τ_hat = params["τ_hat"]
	F_a, τ_a = illuminated_aerodynamics_online(F_hat, τ_hat, ρ_inf, V_inf)

    # Control forces and moments in body_frame
    F_c, τ_c = reaction_control_system(u, params)

	# Forces and moments in body-frame
	F = F_a + F_c
	τ = τ_a + τ_c

	ẋ[1:3] = v # velocity in world frame
	ẋ[4:7] = 0.5 * qmult(q, [0; ω]) #quaternion derivative
	ẋ[8:10] = - g * x / norm(x)  + (1 / m) * qrot(q, F) #acceleration in world frame
	#Euler's equation: I*ω + ω x I*ω = constraint_decrease_ratio
	ẋ[11:13] = Jinv * (τ - cross(ω , J * ω))
	return ẋ
end

function entry_vehicle_simplified_dynamics!(ẋ, X, u, params)
	## States: X ∈ R^13; q = [s;v]
	# x
	# y
	# z
	# q0
	# q1
	# q2
	# q3
	# xdot
	# ydot
	# zdot
	# omega1
	# omega2
	# omega3

	x = X[1:3]
	q = X[4:7]./norm(X[4:7]) #normalize quaternion
	v = X[8:10]
	ω = X[11:13]

	# Parameters
	r_min = params["r_min"] # smallest radius of the cone [m]
	r_max = params["r_max"] # largest radius of the cone [m]
	l = params["l"] # largest radius of the cone [m]
	m = params["m"] # mass of the Phoenix entry system [kg]
	x_g = params["x_g"] # axial center-of-gravity location [m]
	Jxx = params["Jxx"] # Phoenix entry system [kg.m^2]
	Jyy = params["Jyy"] # Phoenix entry system [kg.m^2]
	Jzz = Jyy # Phoenix entry system [kg.m^2]
	Jxy = 0.0 # Phoenix entry system [kg.m^2]
	Jxz = 0.0 # Phoenix entry system [kg.m^2]
	Jyz = 0.0 # Phoenix entry system [kg.m^2]
	g = params["g"] # Mars gravity [m^2.s^-1]
	ρ_0 = params["ρ_0"] # Surface density of Mars atmosphere [kg.m^-3]
	ω_p = [0, 0, 0] # Angular velocity of Mars [rad.s^-1]
	J = Matrix([Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz])# Inertia matrix [kg.m^2]
	Jinv = inv(J) # inverted inertia matrix [kg^-1.m^-2]

	# Mars atmosphere model
	# We simplify the atmosphere model to a constant density.  
	ρ_inf = ρ_0

	# Freestream Velocity
	# Static fluid assumption
	V_inf = - qrot(qconj(q), v + cross(x, ω_p))

	# Compute the force and torque applied of the vehicle at the center of mass.
	F_hat = params["F_hat"]
	τ_hat = params["τ_hat"]
	F_a, τ_a = illuminated_aerodynamics_online(F_hat, τ_hat, ρ_inf, V_inf)

    # Control forces and moments in body_frame
    F_c, τ_c = reaction_control_system(u, params)

	# Forces and moments in body-frame
	F = F_a + F_c
	τ = τ_a + τ_c

	ẋ[1:3] = v # velocity in world frame
	ẋ[4:7] = 0.5 * qmult(q, [0; ω]) #quaternion derivative
	ẋ[8:10] = [g, 0, 0] + (1 / m) * qrot(q, F) #acceleration in world frame
	#Euler's equation: I*ω + ω x I*ω = constraint_decrease_ratio
	ẋ[11:13] = Jinv * (τ - cross(ω , J * ω))
	return ẋ
end

function entry_vehicle_mrp_dynamics!(ẋ, X, u, params)
    ## States: X ∈ R^12; r = modified Rodrigues parameters
    # x
    # y
    # z
    # r1
    # r2
    # r3
    # xdot
    # ydot
    # zdot
    # omega1
    # omega2
    # omega3

    x = X[1:3]
    r = X[4:6] # mrp
    v = X[7:9]
    ω = X[10:12]


	# Parameters
	r_min = params["r_min"] # smallest radius of the cone [m]
	r_max = params["r_max"] # largest radius of the cone [m]
	l = params["l"] # largest radius of the cone [m]
	m = params["m"] # mass of the Phoenix entry system [kg]
	x_g = params["x_g"] # axial center-of-gravity location [m]
	Jxx = params["Jxx"] # Phoenix entry system [kg.m^2]
	Jyy = params["Jyy"] # Phoenix entry system [kg.m^2]
	Jzz = params["Jzz"] # Phoenix entry system [kg.m^2]
	Jxy = params["Jxy"] # Phoenix entry system [kg.m^2]
	Jxz = params["Jxz"] # Phoenix entry system [kg.m^2]
	Jyz = params["Jyz"] # Phoenix entry system [kg.m^2]
	g = params["g"] # Mars gravity [m^2.s^-1]
	h = params["h"] # Atmospheric Scale Height [m]
	ρ_0 = params["ρ_0"] # Surface density of Mars atmosphere [kg.m^-3]
	r_p = params["r_p"] # Volumetric mean radius of Mars [m]
	ω_p = params["ω_p"] # Angular velocity of Mars [rad.s^-1]

	J = Matrix([Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz])# Inertia matrix [kg.m^2]
	Jinv = inv(J) # inverted inertia matrix [kg^-1.m^-2]

    # Rotation
    R = mrptodcm(r)

    # Mars atmosphere model
    # We simplify the atmosphere model to a constant density.  
    ρ_inf = ρ_0

    # Freestream Velocity
    V_inf = - R' * (v + cross(x, ω_p)) 

    # Compute the force and torque applied of the vehicle at the center of mass.
    F_hat = params["F_hat"]
    τ_hat = params["τ_hat"]
    F_a, τ_a = illuminated_aerodynamics_online(F_hat, τ_hat, ρ_inf, V_inf)

    # Control forces and moments in body_frame
    F_c, τ_c = reaction_control_system(u, params)

    # Forces and moments in body-frame
    F = F_a + F_c
    τ = τ_a + τ_c

    ẋ[1:3] = v # velocity in world frame
    ẋ[4:6] = 0.25 * ((1 - r'*r) * ω - 2 * cross(ω,r) + 2 * (ω' * r) * r) #quaternion derivative
    ẋ[7:9] = - g * x / norm(x)  + (1 / m) * qrot(q, F) #acceleration in world frame
    #Euler's equation: I*ω + ω x I*ω = constraint_decrease_ratio
    ẋ[10:12] = Jinv * (τ - cross(ω , J * ω))
    return ẋ
end

function entry_vehicle_mrp_simplified_dynamics!(ẋ, X, u, params)
    ## States: X ∈ R^12; r = modified Rodrigues parameters
    # x
    # y
    # z
    # r1
    # r2
    # r3
    # xdot
    # ydot
    # zdot
    # omega1
    # omega2
    # omega3

    x = X[1:3]
    r = X[4:6] # mrp
    v = X[7:9]
    ω = X[10:12]

    # Parameters
    r_min = params["r_min"] # smallest radius of the cone [m]
    r_max = params["r_max"] # largest radius of the cone [m]
    l = params["l"] # largest radius of the cone [m]
    m = params["m"] # mass of the Phoenix entry system [kg]
    x_g = params["x_g"] # axial center-of-gravity location [m]
    Jxx = params["Jxx"] # Phoenix entry system [kg.m^2]
    Jyy = params["Jyy"] # Phoenix entry system [kg.m^2]
    Jzz = Jyy # Phoenix entry system [kg.m^2]
    Jxy = 0.0 # Phoenix entry system [kg.m^2]
    Jxz = 0.0 # Phoenix entry system [kg.m^2]
    Jyz = 0.0 # Phoenix entry system [kg.m^2]
    g = params["g"] # Mars gravity [m^2.s^-1]
    ρ_0 = params["ρ_0"] # Surface density of Mars atmosphere [kg.m^-3]
    ω_p = [0, 0, 0] # Angular velocity of Mars [rad.s^-1]
    J = Matrix([Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz])# Inertia matrix [kg.m^2]
    Jinv = inv(J) # inverted inertia matrix [kg^-1.m^-2]

    # Rotation
    R = mrptodcm(r)

    # Mars atmosphere model
    # We simplify the atmosphere model to a constant density.  
    ρ_inf = ρ_0

    # Freestream Velocity
    V_inf = - R' * (v + cross(x, ω_p)) 

    # Compute the force and torque applied of the vehicle at the center of mass.
    F_hat = params["F_hat"]
    τ_hat = params["τ_hat"]
    F_a, τ_a = illuminated_aerodynamics_online(F_hat, τ_hat, ρ_inf, V_inf)

    # Control forces and moments in body_frame
    F_c, τ_c = reaction_control_system(u, params)

    # Forces and moments in body-frame
    F = F_a + F_c
    τ = τ_a + τ_c

    ẋ[1:3] = v # velocity in world frame
    ẋ[4:6] = 0.25 * ((1 - r'*r) * ω - 2 * cross(ω,r) + 2 * (ω' * r) * r) #quaternion derivative
    ẋ[7:9] = [g, 0, 0] + (1 / m) * (R * F) #acceleration in world frame

    #Euler's equation: I*ω + ω x I*ω = constraint_decrease_ratio
    ẋ[10:12] = Jinv * (τ - cross(ω , J * ω))
    return ẋ
end

function entry_vehicle_dynamics(X, u, params)
	ẋ = zeros(13,1)
	entry_vehicle_dynamics!(ẋ, X, u, params)
	return ẋ
end

function entry_vehicle_simplified_dynamics(X, u, params)
	ẋ = zeros(13,1)
	entry_vehicle_simplified_dynamics!(ẋ, X, u, params)
	return ẋ
end

function entry_vehicle_mrp_dynamics(X, u, params)
    ẋ = zeros(13,1)
    entry_vehicle_mrp_simplified_dynamics!(ẋ, X, u, params)
    return ẋ
end

function entry_vehicle_mrp_simplified_dynamics(X, u, params)
    ẋ = zeros(13,1)
    entry_vehicle_mrp_simplified_dynamics!(ẋ, X, u, params)
    return ẋ
end


# if PROGRAM_FILE == @__FILE__ 
# 	X_0 = zeros(13, 1)
# 	X_0[1:3] = [3.4e6, 0, 0]
# 	X_0[4:7] = [1, 0, 0, 0]
# 	X_0[8:10] = [-1, 238, 0]
# 	X_0[11:13] = [0, 0, 0]
# 	u_0 = zeros(4, 1)
# 	println("X_0", X_0)
# 	println("u_0", u_0)
# 	println("entry_vehicle_dynamics(X_0,u_0)", entry_vehicle_dynamics(X_0,u_0))
# end

# params = Dict("r_min" => 0.20, # smallest radius of the cone [m]
# 			  "r_max" => 1.30, # largest radius of the cone [m]
#             "δ" => 70 / 360 * 2 * pi, # opening amgle of the cone [rad]
#             "l" => r_max / tan(δ), # length of the cone [m]			  "m" => 569.7, # mass of the Phoenix entry system [kg]
# 			  "x_g" => -0.152, # axial center-of-gravity location [m]
#             "c" => 1.20, # distance between the radial axis and the cluster of thrusters[m]
# 			  "Jxx" => 293.15, # Phoenix entry system [kg.m^2]
# 			  "Jyy" => 184, # Phoenix entry system [kg.m^2]
# 			  "Jzz" => 208.02, # Phoenix entry system [kg.m^2]
# 			  "Jxy" => 0.451, # Phoenix entry system [kg.m^2]
# 			  "Jxz" => -4.424, # Phoenix entry system [kg.m^2] 
# 			  "Jyz" => 0.372, # Phoenix entry system [kg.m^2]
# 			  "g" => 3.711, # Mars gravity [m^2.s^-1]
# 			  "h" => 10800, # Atmospheric Scale Height [m]
# 			  "ρ_0" => 0.020, # Surface density of Mars atmosphere [kg.m^-3]
# 			  "r_p" => 3389.5e3, # Volumetric mean radius of Mars [m]
# 			  "ω_p" => [0, 0, 7.088e-05]) # Angular velocity of Mars [rad.s^-1]
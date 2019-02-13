using LinearAlgebra
using PyCall
using PyPlot
println("FINISHED IMPORTS")

include("AerodynamicsCoefficients.jl")

## Utilities
"""
@(SIGNATURES)
	Conjuguate of a quaternion. 
"""
function qconj(q)
	p = zeros(4)
	p[1] = q[1]
	p[2:4] = - q[2:4]
	return p / norm(q)^2
end

"""
@(SIGNATURES)
	Rotate a vector by a quaternion
"""
function qrot(q,r)
	r + 2*cross(q[2:4],cross(q[2:4],r) + q[1]*r)
end

"""
@(SIGNATURES)
	Multiplication of two quaternions (q = [s;v])
"""
function qmult(q1,q2)
	[q1[1]*q2[1] - q1[2:4]'*q2[2:4]; q1[1]*q2[2:4] + q2[1]*q1[2:4] + cross(q1[2:4],q2[2:4])]
end

function entry_vehicle_dynamics!(ẋ,X,u)
	#TODO change concatenations to make faster!
	# Quaternion representation
	# Modified from D. Mellinger, N. Michael, and V. Kumar,

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
	r_min = 0.20 # smallest radius of the cone [m]
	r_max = 1.30 # largest radius of the cone [m]
	l = 0.55 # length of the cone [m]
	m = 569.7 # mass of the Phoenix entry system [kg]
	x_g = -0.152 # axial center-of-gravity location [m]
	Jxx = 293.15 # Phoenix entry system [kg.m^2]
	Jyy = 184 # Phoenix entry system [kg.m^2]
	Jzz = 208.02 # Phoenix entry system [kg.m^2]
	Jxy = 0.451 # Phoenix entry system [kg.m^2]
	Jxz = -4.424 # Phoenix entry system [kg.m^2]
	Jyz = 0.372 # Phoenix entry system [kg.m^2]
	J = Matrix([Jxx Jxy Jxz; Jxy Jyy Jyz; Jxz Jyz Jzz])# Inertia matrix [kg.m^2]
	Jinv = inv(J) # inverted inertia matrix [kg^-1.m^-2]
	g = 3.711 # Mars gravity [m^2.s^-1]
	h = 10800 # Atmospheric Scale Height [m]
	ρ_0 = 0.020 # Surface density of Mars atmosphere [kg.m^-3]
	r_p = 3389.5e3 # Volumetric mean radius of Mars [m]
	ω_p = [0, 0, 7.088e-05] # Angular velocity of Mars [rad.s^-1]

	# Mars atmosphere model
    altitude = norm(x) - r_p
    ρ_inf = ρ_0 * exp(- altitude / h)

    # Freestream Velocity
    # Static fluid assumption
    V_inf = - qrot(qconj(q), v + cross(x, ω_p))

    # This rotation matrix includes the 1/tan(α) terms required to rescale the aerodynamics 
    # coefficients C_Y2 and C_n2.
    tan_α_cos_θ = V_inf[2] / V_inf[1]
    tan_α_sin_θ = V_inf[3] / V_inf[1]
    Rot_2_1 = [1 0 0; 0 tan_α_cos_θ -tan_α_sin_θ; 0 tan_α_sin_θ tan_α_cos_θ]
    # Since α ∈ [0, π/2[ we can use atan with one variable, it is differentiable in 0.
    ϵ = 1e-100
    α = atan(sqrt(V_inf[2]^2 + V_inf[3]^2 + ϵ^2) - ϵ / V_inf[1])
    
    # Aerodynamic Forces and Moments in Freestream frame
    C_F2, C_τ2 = aerodynamics_coefficents_chebyshev_tan(α, r_max, r_min, l, x_g)

    F_2 = 1.0/2 * ρ_inf * norm(V_inf)^2 * C_F2
    τ_2 = 1.0/2 * ρ_inf * norm(V_inf)^2 * C_τ2

    F_a = Rot_2_1 * F_2
    τ_a = Rot_2_1 * τ_2

    # Control forces and moments in body_frame
    F_c = zeros(3)
    τ_c = zeros(3)

    # Forces and moments in body-frame
    F = F_a + F_c
    τ = τ_a + τ_c

    ẋ[1:3] = v # velocity in world frame
    ẋ[4:7] = 0.5 * qmult(q, [0; ω]) #quaternion derivative
    ẋ[8:10] = - g * x / norm(x)  + (1 / m) * qrot(q, F) #acceleration in world frame
	#Euler's equation: I*ω + ω x I*ω = constraint_decrease_ratio
    ẋ[11:13] = Jinv * (τ - cross(ω , J * ω))
    return ẋ, V_inf
end

function entry_vehicle_dynamics(X,u)
	ẋ = zeros(13,1)
	entry_vehicle_dynamics!(ẋ,X,u)
	ẋ
end

if PROGRAM_FILE == @__FILE__ 
	X_0 = zeros(13, 1)
	X_0[1:3] = [3.4e6, 0, 0]
	X_0[4:7] = [1, 0, 0, 0]
	X_0[8:10] = [-1, 238, 0]
	X_0[11:13] = [0, 0, 0]
	u_0 = zeros(4, 1)
	println("X_0", X_0)
	println("u_0", u_0)
	println("entry_vehicle_dynamics(X_0,u_0)", entry_vehicle_dynamics(X_0,u_0))
end

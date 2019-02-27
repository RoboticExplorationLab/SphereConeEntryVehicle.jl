using HCubature
using LinearAlgebra
using StaticArrays

function cross_product_matrix(a)
    a_mat = [0 -a[3] a[2]; a[3] 0 -a[1]; -a[2] a[1] 0]
    return a_mat
end

function illuminated_aerodynamics(r_max, r_min, l, x_g, ρ_inf, V_inf)
    F_hat, τ_hat = illuminated_integral_aerodynamics(r_max, r_min, l, x_g, ρ_inf, V_inf)
    F_a = ρ_inf * F_hat * kron(V_inf, V_inf)
    τ_a = ρ_inf * τ_hat * kron(V_inf, V_inf)
    return F_a, τ_a
end

function illuminated_aerodynamics_online(F_hat, τ_hat, ρ_inf, V_inf)
    F_a = ρ_inf * F_hat * kron(V_inf, V_inf)
    τ_a = ρ_inf * τ_hat * kron(V_inf, V_inf)
    return F_a, τ_a
end

function illuminated_integral_aerodynamics(r_max, r_min, l, x_g)
    F_hat_cone, τ_hat_cone = illuminated_integral_aerodynamics_cone(r_max, r_min, l)
    F_hat_sphere, τ_hat_sphere = illuminated_integral_aerodynamics_sphere(r_max, r_min, l)
    
    # Geometry of the vehicle
    # angle of the cone
    δ = atan(r_max, l)
    # radius of the sphere
    r_d = r_min / cos(δ)
    # distance between the center of the sphere and the base of the small cone
    r_0 = sqrt(r_d^2 - r_min^2)
    # length of the section of cone.
    f = l * (1 - r_min/r_max)
    GC = [0, 0, 0] - [x_g, 0, 0]
    GS = [f - r_0, 0, 0] - [x_g, 0, 0]
    GC_mat = cross_product_matrix(GC)
    GS_mat = cross_product_matrix(GS)   
    F_hat = F_hat_cone + F_hat_sphere
    τ_hat = τ_hat_cone + GC_mat * F_hat_cone + τ_hat_sphere + GS_mat * F_hat_sphere
    return F_hat, τ_hat
end

function illuminated_aerodynamics_cone_only(r_max, r_min, l, x_g, ρ_inf, V_inf)
    F_hat, τ_hat = illuminated_integral_aerodynamics_cone_only(r_max, r_min, l, x_g)
    F_a = ρ_inf * F_hat * kron(V_inf, V_inf)
    τ_a = ρ_inf * τ_hat * kron(V_inf, V_inf)
    return F_a, τ_a
end

function illuminated_integral_aerodynamics_cone_only(r_max, r_min, l, x_g)
    F_hat_cone, τ_hat_cone = illuminated_integral_aerodynamics_cone(r_max, r_min, l)

    # Geometry of the vehicle
    # angle of the cone
    δ = atan(r_max, l)
    # radius of the sphere
    # r_d = r_min / cos(δ)
    # distance between the center of the sphere and the base of the small cone
    # r_0 = sqrt(r_d^2 - r_min^2)
    # length of the section of cone.
    f = l * (1 - r_min/r_max)
    GC = [0, 0, 0] - [x_g, 0, 0]
    # GS = [f - r_0, 0, 0] - [x_g, 0, 0]
    GC_mat = cross_product_matrix(GC)
    # GS_mat = cross_product_matrix(GS)   
    F_hat = F_hat_cone #+ F_hat_sphere
    τ_hat = τ_hat_cone + GC_mat * F_hat_cone #+ τ_hat_sphere + GS_mat * F_hat_sphere
    return F_hat, τ_hat
end

function cone_force_integrand_parameterized(x::SVector{2,Float64}, r_max, r_min, l)
    u = x[1]
    v = x[2]
    # angle of the cone
    δ = atan(r_max, l)
    du = [- 1/tan(δ), cos(v), -sin(v)]
    dv = [0, -u*sin(v), -u*cos(v)]
    ds = norm(cross(du, dv))
    n = cross(du, dv) / ds
    y = n * vec(n * n')' * ds
    return y
end

function cone_moment_integrand_parameterized(x::SVector{2,Float64}, r_max, r_min, l)
    u = x[1]
    v = x[2]
    # angle of the cone
    δ = atan(r_max, l)
    du = [- 1/tan(δ), cos(v), -sin(v)]
    dv = [0, -u*sin(v), -u*cos(v)]
    r = [l - u/tan(δ), u*cos(v), -u*sin(v)]
    ds = norm(cross(du, dv))
    n = cross(du, dv) / ds
    y = cross(r, n) * vec(n * n')' * ds
    return y
end

function sphere_force_integrand_parameterized(x::SVector{2,Float64}, r_d)
    u = x[1]
    v = x[2]
    du = [1, -u*cos(v)/sqrt(r_d^2 - u^2), -u*sin(v)/sqrt(r_d^2 - u^2)]
    dv = [0, -sqrt(r_d^2 - u^2)*sin(v), sqrt(r_d^2 - u^2)*cos(v)]
    ds = norm(cross(du, dv))
    n = cross(du, dv) / ds
    y = n * vec(n * n')' * ds
    return y
end

function sphere_moment_integrand_parameterized(x::SVector{2,Float64}, r_d)
    u = x[1]
    v = x[2]
    du = [1, -u*cos(v)/sqrt(r_d^2 - u^2), -u*sin(v)/sqrt(r_d^2 - u^2)]
    dv = [0, -sqrt(r_d^2 - u^2)*sin(v), sqrt(r_d^2 - u^2)*cos(v)]
    r = [u, sqrt(r_d^2 - u^2)*cos(v), sqrt(r_d^2 - u^2)*sin(v)]
    ds = norm(cross(du, dv))
    n = cross(du, dv) / ds
    y = cross(r, n) * vec(n * n')' * ds
    return y
end

function illuminated_integral_aerodynamics_cone(r_max, r_min, l)
    # define domain of integration and number of evaluation in the integral computation.
    start_bound = [r_min, 0.0]
    end_bound = [r_max, 2 * pi]
    # Define both integrands with fixed parameters r_d and α.
    function cone_force_integrand(x::SVector{2,Float64})
        return cone_force_integrand_parameterized(x, r_max, r_min, l)
    end
    function cone_moment_integrand(x::SVector{2,Float64})
        return cone_moment_integrand_parameterized(x, r_max, r_min, l)
    end
    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000
    # Integration
    F_hat, force_error = hcubature(cone_force_integrand, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    τ_hat, moment_error = hcubature(cone_moment_integrand, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    return F_hat, τ_hat
end

function illuminated_integral_aerodynamics_sphere(r_max, r_min, l)
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
    function sphere_force_integrand(x::SVector{2,Float64})
        return sphere_force_integrand_parameterized(x, r_d)
    end
    function sphere_moment_integrand(x::SVector{2,Float64})
        return sphere_moment_integrand_parameterized(x, r_d)
    end
    # Parameters of integration to limit memory usage and running time.
    rtol = 1e-4
    maxevals = 1000
    # Integration
    F_hat, force_error = hcubature(sphere_force_integrand, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    τ_hat, moment_error = hcubature(sphere_moment_integrand, start_bound, end_bound;
        norm=norm, rtol=rtol, maxevals=maxevals)
    return F_hat, τ_hat
end


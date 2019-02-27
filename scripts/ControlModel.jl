using LinearAlgebra

function reaction_control_system(u, params)
	# u = [T11, T12, T13, T21, T22, T23, T31, T32, T33, T41, T42, T43]
	T11, T12, T13, T21, T22, T23, T31, T32, T33, T41, T42, T43 = u
	# Forces for the 4 clusrers of thrusters.
	F = [T11+T12  0     -T13;
		 T21+T22  -T23  0   ;
		 T31+T32  0     T33 ;
		 T41+T42  T43   0   ] 

	x_g = params["x_g"]
	c = params["c"]
	# Positions of the clusters 
	C = [0 +c 0 ; # C_1
		 0 0  +c; # C_2
		 0 -c 0 ; # C_3
		 0 0  -c] # C_4
	# Position of the center of mass.            
	G = [x_g, 0, 0]

	F_c = reshape(sum(F, dims=1), 3)
	τ_c = zeros(3)
	for i=1:4
		GC_i = C[i, :] - G
		τ_c += cross(GC_i, F[i, :])
	end   
	return F_c, τ_c
end


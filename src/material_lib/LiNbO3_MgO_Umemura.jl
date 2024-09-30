export MgO_LiNbO₃_Umemura

"""
These functions create a symbolic representation (using ModelingToolkit) of the
Sellmeier Equation model for the temperature and wavelength dependence
of 5% MgO:LiNbO3's ordinary and extraordinary indices of refraction for 0.4133μm < λ < 2.03μm
Note that these functions are for MgO-doped congruent LiNbO3 (cLN), not stoichiometric LiNbO3 (sLN).
Equation form is based on "Temperature-dependent phase-matching properties with oo-e and oo-o 
interactions in 5 mol.% MgO-doped congruent LiNbO3" by Matsuda, ... Umemura et al. (2015)
https://doi.org/10.1117/12.2078801

This model is then exported to other functions that use it and its
derivatives to return index, group index and GVD values as a function
of temperature and wavelength.
Variable units are λ in [um] and T in [ᵒC]
"""

# # Thermo-optic coefficients
# function dnₒ_dT(λ, T)
#     (0.5017/λ^2 - 0.4290/λ + 0.3954) * 10.0^(-5) * (1 + 0.00216*T)
# end

# function dnₑ_dT(λ, T)
#     (0.4175/λ^3 - 0.6643/λ^2 + 0.9036/λ + 3.5332 - 0.0744*λ) * 10.0^(-5) * (1 + 0.00276*T)
# end

# # Integrated form of above equations
# function Δnₒ(λ, T)
#     ΔT = T-T₀
#     (0.5017/λ^2 - 0.4290/λ + 0.3954) * 10.0^(-5) * (ΔT + 0.00108*(ΔT)^2)
# end

# function Δnₑ(λ, T)
#     ΔT = T-T₀
#     (0.4175/λ^3 - 0.6643/λ^2 + 0.9036/λ + 3.53320 - 0.0744λ) * 10.0^(-5) * (ΔT + 0.00138*(ΔT)^2)
# end

# function nₒ²_sym(λ, T)
#     nₒ_at_T₀ = (4.8898 + 0.11738/(λ^2 - 0.04564) - 0.02662*λ^2)^0.5
#     (nₒ_at_T₀ + Δnₒ(λ, T))^2
# end

# function nₑ²_sym(λ, T)
#     nₑ_at_T₀ = (4.54514 + 0.096471/(λ^2 - 0.043763) - 0.021502*λ^2)^0.5
#     (nₑ_at_T₀ + Δnₑ(λ, T))^2
# end

# TODO - hopefully function subcall plays nice with symbolic.jl upon use of substitute.

T₀ = 20.0 # reference temperature in [Deg C]

function dnₒ_dT_ω(ω, T)
    (0.5017*ω^2 - 0.4290*ω + 0.3954) * 10.0^(-5) * (1 + 0.00216*T)
end

function dnₑ_dT_ω(ω, T)
    (0.4175*ω^3 - 0.6643*ω^2 + 0.9036*ω + 3.5332 - 0.0744/ω) * 10.0^(-5) * (1 + 0.00276*T)
end

function Δnₒ_ω(ω, T)
    ΔT = T-T₀
    (0.5017*ω^2 - 0.4290*ω + 0.3954) * 10.0^(-5) * (ΔT + 0.00108*(ΔT)^2)
end

function Δnₑ_ω(ω, T)
    ΔT = T-T₀
    (0.4175*ω^3 - 0.6643*ω^2 + 0.9036*ω + 3.53320 - 0.0744/ω) * 10.0^(-5) * (ΔT + 0.00138*(ΔT)^2)
end

function nₒ²_sym_ω(ω, T)
    nₒ_at_T₀ = (4.8898 + (0.11738*ω^2) /(1 - 0.04564*ω^2) - 0.02662/ω^2)^(1/2)
    (nₒ_at_T₀ + Δnₒ_ω(ω, T))^2
end

function nₑ²_sym_ω(ω, T)
    nₑ_at_T₀ = (4.54514 + (0.096471*ω^2)/(1 - 0.043763*ω^2) - 0.021502/ω^2)^(1/2)
    (nₑ_at_T₀ + Δnₑ_ω(ω, T))^2
end

pᵪ₂ = (
	d₃₃ =   20.3,    #   pm/V
	d₃₁ =   -4.1,    #   pm/V
	d₂₂ =   2.1,     #   pm/V
	λs  =  [1.313, 1.313, 1.313/2.0]
)

function make_MgO_LiNbO₃(;pₒ=pₒ,pₑ=pₑ,pᵪ₂=pᵪ₂)
	@variables λ, ω, T, λs[1:3]
	# nₒ² = n²_MgO_LiNbO₃_sym(λ, T; pₒ...)
	# nₑ² = n²_MgO_LiNbO₃_sym(λ, T; pₑ...)
	nₒ² = nₒ²_sym_ω(ω, T)
	nₑ² = nₑ²_sym_ω(ω, T)
	ε 	= diagm([nₒ², nₒ², nₑ²])
	d₃₃, d₃₁, d₂₂, λᵣs = pᵪ₂
	χ⁽²⁾ᵣ = cat(
		[ 	0.0	 	-d₂₂ 	d₃₁			#	xxx, xxy and xxz
		 	-d₂₂	0.0 	0.0			#	xyx, xyy and xyz
			d₃₁	 	0.0		0.0		],	#	xzx, xzy and xzz
		[ 	-d₂₂	0.0 	0.0			#	yxx, yxy and yxz
			0.0	 	d₂₂ 	d₃₁			#	yyx, yyy and yyz
			0.0	 	d₃₁		0.0		],	#	yzx, yzy and yzz
		[ 	d₃₁	 	0.0 	0.0			#	zxx, zxy and zxz
			0.0	 	d₃₁ 	0.0			#	zyx, zyy and zyz
			0.0	 	0.0 	d₃₃		],	#	zzx, zzy and zzz
		 dims = 3
	)
	
	# ngₒ = ng_model(nₒ,ω)
	# gvdₒ = gvd_model(nₒ,ω)
	
	ε_λ = substitute.(ε,(Dict([(ω=>1/λ),]),))
	nₒ_λ,nₑ_λ = sqrt.(substitute.((nₒ²,nₑ²),(Dict([(ω=>1/λ),]),)))
	# nₑ = sqrt(nₑ²)
	# nₒ = sqrt(nₒ²)
	# ngₑ = ng_model(nₑ,ω)
	# gvdₑ = gvd_model(nₑ,ω)
	models = Dict([
		:nₒ		=>	nₒ_λ,
		# :ngₒ	=>	ngₒ,
		# :gvdₒ	=>	gvdₒ,
		:nₑ		=>	nₑ_λ,
		# :ngₑ	=>	ngₑ,
		# :gvdₑ	=>	gvdₑ,
		# :ng		=>	diagm([ngₒ, ngₒ, ngₑ]),
		# :gvd	=>	diagm([gvdₒ, gvdₒ, gvdₑ]),
		:ε 		=> 	ε,
		:χ⁽²⁾	=>	SArray{Tuple{3,3,3}}(Δₘ(λs,ε_λ, λᵣs, χ⁽²⁾ᵣ)),
	])
	defaults =	Dict([
		:ω		=>		inv(0.8),	# μm⁻¹
		:λ		=>		0.8,		# μm
		:T		=>		24.5,		# °C
		:λs₁	=>		1.064,		# μm
		:λs₂	=>		1.064,		# μm
		:λs₃	=>		0.532,		# μm

	])
	Material(
		models,
		defaults,
		:MgO_LiNbO₃_Umemura,
		colorant"seagreen2", #RGB{N0f8}(0.22,0.596,0.149),	#
		)
end

################################################################

MgO_LiNbO₃_Umemura = make_MgO_LiNbO₃()


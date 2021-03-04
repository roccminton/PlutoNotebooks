### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 728d5546-7a96-11eb-1757-93681c4e7955
begin
	import Pkg
	
	Pkg.add("PlutoUI")
	Pkg.add("Plots")
	
	using PlutoUI
	using Plots
end

# ╔═╡ 58ba05d0-7a96-11eb-08da-41b782596598
md"# Euler Method"

# ╔═╡ 97f8a9f8-7b78-11eb-2f25-e5228f41db56
md"We implement a Euer Method to solve first oderd ODEs and use it to some relevant Biological Models."

# ╔═╡ ecad4e78-7a96-11eb-2b54-e7f910e83835
function euler(f,t0,tn,Δt,x₀)
	T = t0:Δt:tn
	F = Vector{typeof(x₀)}(undef,length(T))
	F[1] = x₀
	for (n,t) ∈ enumerate(T[2:end])
		F[n+1] = F[n] .+ f(t,F[n]) .* Δt
	end
	return F
end

# ╔═╡ cd0f0650-7b78-11eb-12c0-cb427957850b
md"""
---
### Logistic Population Growth
"""

# ╔═╡ e9b52262-7b7d-11eb-21aa-2bf0d9f5ef50
md"""

$\frac{dP}{dt} = rP \left( 1- \frac{P}{K} \right)$

"""

# ╔═╡ 04fdd000-7b79-11eb-029b-b3765364f627
md"""
growth rate r $(@bind r NumberField(0.0001:0.0001:0.1,default=0.01))

carrying capacity K $(@bind K NumberField(1:10^7,default=10_000))
"""

# ╔═╡ dcffb046-7b78-11eb-333c-7df614c12617
begin
	LogisticGrowth(t,x) = r*x*(1-x/K)
	LG_0 = 1.0
end;

# ╔═╡ 5b3c4954-7b79-11eb-11ba-57b8ae51a6f5
let
	Δt = 0.01
	
	p = plot(framestlye=:zerolines,legend=:topleft)
	
	T=0:Δt:2500
	F = euler(LogisticGrowth,T[1],T[end],Δt,LG_0)
	plot!(p,T,F,label="Population Size")
	
	p
end

# ╔═╡ 7a010972-7b78-11eb-3478-f1fbc1f0645d
md"""
---
### Lotka-Volterra
"""

# ╔═╡ b963fdb2-7b7e-11eb-3184-1b216f77bd64
md"""

${\frac{dN_1}{dt} = N_1 (\varepsilon_1 - \gamma_1 N_2) \quad , \quad \frac{dN_2}{dt} = -N_2 (\varepsilon_2 - \gamma_2 N_1)}$ 


"""

# ╔═╡ 07b392ee-7b77-11eb-00be-876772164be8
md"""
Reproduktionsraten
	ε₁ $(@bind ε₁ Slider(0.01:0.01:1,show_value=true, default=0.5))
	ε₂ $(@bind ε₂ Slider(0.01:0.01:1,show_value=true, default=0.25))

Sterberaten 
	δ₁ $(@bind δ₁ Slider(0.01:0.01:1,show_value=true, default=0.07))
	δ₂ $(@bind δ₂ Slider(0.01:0.01:1,show_value=true, default=0.07))
"""

# ╔═╡ 7fc751be-7a96-11eb-207d-d5d0142f43ba
begin
	LotkaVolterra(t,x) = [
			x[1]*(ε₁-δ₁*x[2]),
			-x[2]*(ε₂-δ₂*x[1])
	]
	LV_0 = [10.0,10.0]
end;

# ╔═╡ 62645b8c-7b7a-11eb-1b35-a3957a0252fb
md"""
---
### SIR-Model
"""

# ╔═╡ 4341adf4-7b7f-11eb-25c1-e546febcad3e
md"""

${\frac{dS}{dt} = \nu N - \beta \frac{SI}{N} - \mu S}$

${\frac{dI}{dt} = \beta \frac{SI}{N} - \gamma I - \mu I}$

${\frac{dR}{dt} = \gamma I - \mu R \phantom{-\beta \frac{SI}{N}}}$

"""

# ╔═╡ e64e9d62-7b7b-11eb-3620-a7a57b954fcc
md"""
recovery rate γ $(@bind γ Slider(0.0:0.01:1,show_value=true,default=0.15)) | 
infection rate β $(@bind β Slider(0.0:0.01:1,show_value=true,default=0.7))

death rate μ $(@bind μ Slider(0.0:0.01:1,show_value=true)) |
birth rate ν $(@bind ν Slider(0.0:0.01:1,show_value=true))

Reduce infection rate if incidence is high by factor of: $(@bind high_inc_f NumberField(0:0.1:100,default=2)) | 
High incidence at $(@bind high_inc NumberField(0:1000,default=50)) per 10^5

"""

# ╔═╡ 6b37eb20-7b7a-11eb-07de-a1042efa6721
begin
	function SIR(t,x)
		N = sum(x)
		S, I, R = x[1], x[2], x[3]
		β = Β(I,N)
		return [
			(ν*N - β*S*I/N - μ*S),	# dS/dt
			(β*S*I/N - γ*I - μ*I),	# dI/dt
			(γ*I - μ*R)				# dR/dt
		]
	end
	
	Β(I,N) = I/N <= high_inc/100_000 ? β/high_inc_f : β
	
	SIR_0 = [9999.9,0.1,0.0]
	
	Δt = 0.01
	T=0:Δt:150
end;

# ╔═╡ 79a2aad2-7c2a-11eb-00b0-6bebd6658619
md"""
---
### Six dimensional SIR-Model
"""

# ╔═╡ 968af7ba-7c2a-11eb-0aa4-6d51510020b1
md"""
Devide the population in two subgroups with different infection rates within and between these groups.

${\frac{dS_1}{dt} = - \beta_{11} \frac{S_1I_1}{N} - \beta_{12} \frac{S_1I_2}{N} \phantom{-\beta \frac{SI}{N}}}$

${\frac{dI_1}{dt} = \phantom{-} \beta_{11} \frac{S_1I_1}{N} + \beta_{12} \frac{S_1I_2}{N} -   \gamma_1 I_1 }$

${\frac{dR_1}{dt} = \phantom{- \beta_{11} \frac{S_1I_1}{N} - \beta_{12} \frac{S_1I_2}{N} -}\gamma_1 I_1 }$

${\frac{dS_2}{dt} = - \beta_{22} \frac{S_2I_2}{N} - \beta_{21} \frac{S_2I_1}{N} \phantom{-\beta \frac{SI}{N}}}$

${\frac{dI_2}{dt} = \phantom{-} \beta_{22} \frac{S_2I_2}{N} + \beta_{21} \frac{S_2I_1}{N} -   \gamma_2 I_2 }$

${\frac{dR_2}{dt} = \phantom{- \beta_{11} \frac{S_1I_1}{N} - \beta_{12} \frac{S_1I_2}{N} -}\gamma_2 I_2 }$


Here $\beta_{11}$ is the infection rate within group one, $\beta_{12}$ is the infection rate from group one to group two, and so on.
"""

# ╔═╡ cbf8bba4-7c26-11eb-3301-07233740f82a
md"""

β₁₁ = $(@bind β11 Slider(0.0:0.01:1.0, show_value=true, default=0.75)) |
β₁₂ = $(@bind β12 Slider(0.0:0.01:1.0, show_value=true, default=0.5)) 

β₂₁ = $(@bind β21 Slider(0.0:0.01:1.0, show_value=true, default=0.5)) |
β₂₂ = $(@bind β22 Slider(0.0:0.01:1.0, show_value=true, default=0.75)) 

γ₁ = $(@bind γ1 Slider(0.0:0.01:1.0, show_value=true, default=0.25)) |
γ₂ = $(@bind γ2 Slider(0.0:0.01:1.0, show_value=true, default=0.25)) 

"""

# ╔═╡ 58e3898c-7c21-11eb-3364-9168fcedb22c
begin
	function SIR6D(t,x)
		N = sum(x)
		S1, S2 = x[1], x[4]
		I1, I2 = x[2], x[5]
		R1, R2 = x[3], x[6]
		
		return [
			(-β11*S1*I1/N - β21*S1*I2/N),			# dS₁/dt
			(β11*S1*I1/N + β21*S1*I2/N - γ1*I1),	# dI₁/dt
			(γ1*I1),								# dR₁/dt
			(-β22*S2*I2/N - β12*S2*I1/N),			# dS₂/dt
			(β22*S2*I2/N + β12*S2*I1/N - γ2*I2),	# dI₂/dt
			(γ2*I2),								# dR₂/dt
		]
	end
	
	SIR6D_0 = [
		4999.0,1.0,0.0,
		5000.0,0.0,0.0,
	]
end;

# ╔═╡ a7ec0730-7b76-11eb-31ec-fdb09304eb10
reshape_result(F) = [[f[i] for f ∈ F] for i ∈ 1:length(F[1])] ;

# ╔═╡ c080dd1a-7a96-11eb-345c-e364dc0f8ac9
let
	Δt = 0.01
	
	p = plot(framestlye=:zerolines)
	
	T=0:Δt:100
	F = euler(LotkaVolterra,T[1],T[end],Δt,LV_0)
	plot!(p,T,reshape_result(F),label=["N₁" "N₂"])
	
	p
end

# ╔═╡ fa94e89e-7c1e-11eb-08a7-91210bc2fe60
begin
	F = euler(SIR,T[1],T[end],Δt,SIR_0)
	rF = reshape_result(F)
	N = sum.(F)
end;

# ╔═╡ 68f66db2-7b7c-11eb-138d-2f40b3bfbb4e
let
	p = plot(framestlye=:zerolines,size=(680,450))
	
	plot!(p,T,rF,label=["S" "I" "R"])
	
	high_inc_f ≠ 1 && plot!(p,T,N .* high_inc/100_000, label="",linestyle=:dash) 
	
	p
end

# ╔═╡ 0ddf4fe6-7c28-11eb-2969-7b856ff408be
begin
	F6D = euler(SIR6D,T[1],T[end],Δt,SIR6D_0)
	rF6D = reshape_result(F6D)
	N6D = sum.(F6D)
end;

# ╔═╡ d1644026-7c27-11eb-17ae-67c414382cdc
let
	p = plot(framestlye=:zerolines,size=(680,450))
	
	plot!(p,T,rF6D,
		#label=["S₁" "I₁" "R₁" "S₂" "I₂" "R₂"],
		label="",
		alpha=0.5,
		lw=1,
		color=[:darkblue :darkred :darkgreen :lightblue :orangered :lightgreen]
	)
	
	plot!(p,T,[rF6D[i] .+ rF6D[i+3] for i ∈ 1:3],
		label=["S" "I" "R"],
		color=[:blue :red :green],
		lw=2
	)
	
	p
end

# ╔═╡ Cell order:
# ╟─58ba05d0-7a96-11eb-08da-41b782596598
# ╟─97f8a9f8-7b78-11eb-2f25-e5228f41db56
# ╠═ecad4e78-7a96-11eb-2b54-e7f910e83835
# ╟─cd0f0650-7b78-11eb-12c0-cb427957850b
# ╟─e9b52262-7b7d-11eb-21aa-2bf0d9f5ef50
# ╠═dcffb046-7b78-11eb-333c-7df614c12617
# ╟─04fdd000-7b79-11eb-029b-b3765364f627
# ╟─5b3c4954-7b79-11eb-11ba-57b8ae51a6f5
# ╟─7a010972-7b78-11eb-3478-f1fbc1f0645d
# ╟─b963fdb2-7b7e-11eb-3184-1b216f77bd64
# ╠═7fc751be-7a96-11eb-207d-d5d0142f43ba
# ╟─07b392ee-7b77-11eb-00be-876772164be8
# ╟─c080dd1a-7a96-11eb-345c-e364dc0f8ac9
# ╟─62645b8c-7b7a-11eb-1b35-a3957a0252fb
# ╟─4341adf4-7b7f-11eb-25c1-e546febcad3e
# ╠═6b37eb20-7b7a-11eb-07de-a1042efa6721
# ╟─fa94e89e-7c1e-11eb-08a7-91210bc2fe60
# ╟─e64e9d62-7b7b-11eb-3620-a7a57b954fcc
# ╟─68f66db2-7b7c-11eb-138d-2f40b3bfbb4e
# ╟─79a2aad2-7c2a-11eb-00b0-6bebd6658619
# ╟─968af7ba-7c2a-11eb-0aa4-6d51510020b1
# ╠═58e3898c-7c21-11eb-3364-9168fcedb22c
# ╟─cbf8bba4-7c26-11eb-3301-07233740f82a
# ╟─d1644026-7c27-11eb-17ae-67c414382cdc
# ╟─0ddf4fe6-7c28-11eb-2969-7b856ff408be
# ╟─728d5546-7a96-11eb-1757-93681c4e7955
# ╟─a7ec0730-7b76-11eb-31ec-fdb09304eb10

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
    Pkg.activate(".")
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
begin
	function euler(f,t0,tn,Δt,x₀)
		T = t0:Δt:tn
		F = Vector{typeof(x₀)}(undef,length(T))
		F[1] = x₀
		for (n,t) ∈ enumerate(T[2:end])
			F[n+1] = F[n] .+ f(t,F[n]) .* Δt
		end
		return F
	end
	
	reshape_result(F) = [[f[i] for f ∈ F] for i ∈ 1:length(F[1])] ;
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

# ╔═╡ c080dd1a-7a96-11eb-345c-e364dc0f8ac9
let
	Δt = 0.01
	
	p = plot(framestlye=:zerolines)
	
	T=0:Δt:100
	F = euler(LotkaVolterra,T[1],T[end],Δt,LV_0)
	plot!(p,T,reshape_result(F),label=["N₁" "N₂"])
	
	p
end

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
"""

# ╔═╡ 6b37eb20-7b7a-11eb-07de-a1042efa6721
begin
	function SIR(t,x)
		N = sum(x)
		S, I, R = x[1], x[2], x[3]
		return [
			(ν*N - β*S*I/N - μ*S),	# dS/dt
			(β*S*I/N - γ*I - μ*I),	# dI/dt
			(γ*I - μ*R)				# dR/dt
		]
	end
	
	SIR_0 = [9999.9,0.1,0.0]
	
	Δt = 0.01
	T=0:Δt:150
end;

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
# ╟─728d5546-7a96-11eb-1757-93681c4e7955

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

# ╔═╡ 8dfd7e3e-80b8-11eb-0a8b-99c7a0102380
md"""# SIR-Models

A collection of variations of SIR Models:

* Standard SIR with birth & death
* Inzidenz dependent infection rates
* 6D-SIR with two subpopulations

First set the time horizon and the step size for the euler method for all simulations below.

Start $(@bind T0 NumberField(0:10_000, default=0)) |  
Step size Δt $(@bind Δt NumberField(0.0001:0.0001:1.0,default=0.01)) |
Stopp $(@bind Tend NumberField(0:10_000, default=150))

"""

# ╔═╡ 62645b8c-7b7a-11eb-1b35-a3957a0252fb
md"""
---
### Classical SIR-Model 

The standard SIR-Model with constant birth and death, infection and recovery rates.


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

# ╔═╡ 81608c56-80b9-11eb-00f0-3b50905ec998
md"""
---
### Inzidenz dependent infection rate

Here the infection rate $\beta=\beta(I)$ is dependent on the infection rate in a population of constant size.

${\frac{dS}{dt} = - \beta(I) \frac{SI}{N} \phantom{- \gamma I} }$

${\frac{dI}{dt} = \phantom{-} \beta(I) \frac{SI}{N} - \gamma I}$

${\frac{dR}{dt} = \phantom{-\beta(I) \frac{SI}{N}-} \gamma I }$

with the following infection rate functions

* Constant: $\beta(I)\equiv\beta$
* Step Function: $\beta(I) = \beta (1 - \mathbb{1}_{[I_{c},\infty)}(I) (1-\frac{1}{r_s}))$
* Logistic Decay: $\beta(I) = \frac{\beta}{1+\exp(r_l(I-I_c))}$

where $I_c$ is the critical incidence at which the measurements stat and $r_s, r_f$ is the factor resp. rate at which the infection rate decreases.
"""

# ╔═╡ 5e81dd2c-80b9-11eb-364e-63fa4ebd27da
md"""
Choice of Incidence Function $(@bind Β_name Select(["LogisticDecay", "StepFunction","Constant"]))
"""

# ╔═╡ 4f8dee16-80ba-11eb-1d0f-f3ecdc9ac671
if Β_name == "LogisticDecay"
	md"""
	maximum infection rate β $(@bind b Slider(0.0:0.01:1,show_value=true,default=0.7))
	
	Maximum decrease of infection rate at $(@bind high_inc NumberField(0:1000,default=50)) per 10^5 | Rate of decrease $(@bind r_ld Slider(0.0:0.001:0.1,show_value=true,default=0.01))
	"""
elseif Β_name == "StepFunction"
	md"""
	high infection rate β $(@bind b Slider(0.0:0.01:1,show_value=true,default=0.7))
	
	Reduce infection rate if incidence is high by factor of: $(@bind high_inc_f NumberField(0:0.1:100,default=2)) | 
	High incidence at $(@bind high_inc NumberField(0:1000,default=50)) per 10^5
	"""
elseif Β_name == "Constant"
	md"""
	infection rate β $(@bind b Slider(0.0:0.01:1,show_value=true,default=0.7))
	"""
else
	md""
end

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

# ╔═╡ 32f4bc7c-80b9-11eb-0111-074cc68e892f
md"---"

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
	
	reshape_result(F) = [[f[i] for f ∈ F] for i ∈ 1:length(F[1])]
	
	T = T0:Δt:Tend
end;

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

# ╔═╡ 10a37328-80bb-11eb-1899-c9abbcfaaaf0
begin
	function SIR_β(t,x)
		S, I, R = x[1], x[2], x[3]
		β = Β(I,N_β)
		return [
			(- β*S*I/N_β ),		# dS/dt
			(β*S*I/N_β - γ*I),	# dI/dt
			(γ*I)				# dR/dt
		]
	end
	
	#Infection Rate Functions
	StepFunction(I,N) = I/N >= high_inc/100_000 ? b/high_inc_f : b
	LogisticDecay(I,N) = b / (1+exp(r_ld*(((I/N)*10^5-high_inc))))
	Constant(I,N) = b
	
	function_choice = Dict(
		"StepFunction" => StepFunction,
		"LogisticDecay" => LogisticDecay,
		"Constant" => Constant
	)
	
	#---
	
	Β = function_choice[Β_name]
	
	SIR_β_0 = [9999.9,0.1,0.0]
	N_β = sum(SIR_β_0)
	
	F_β = euler(SIR_β,T[1],T[end],Δt,SIR_β_0)
	rF_β = reshape_result(F_β)
end;

# ╔═╡ 02299d2a-8010-11eb-024a-650734491c82
let
	p = plot(framestlye=:zerolines,size=(680,450))
	
	plot!(p,T,(rF_β[2]./N_β) .* 10^5,label="Incidence per 100.000",color=:red,legend=:topright) 
	p2 = twinx()
	plot!(p2,T,Β.(rF_β[2],N_β), label="Infection rate",color=:orange,legend=:bottomright)
	
	p
end

# ╔═╡ 90a78c38-80c2-11eb-3124-1dadb02c2806
let
	p = plot(framestlye=:zerolines,size=(680,450))
	
	plot!(p,T,rF_β,label=["S" "I" "R"])
	
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
# ╟─8dfd7e3e-80b8-11eb-0a8b-99c7a0102380
# ╟─62645b8c-7b7a-11eb-1b35-a3957a0252fb
# ╠═6b37eb20-7b7a-11eb-07de-a1042efa6721
# ╟─e64e9d62-7b7b-11eb-3620-a7a57b954fcc
# ╠═68f66db2-7b7c-11eb-138d-2f40b3bfbb4e
# ╟─81608c56-80b9-11eb-00f0-3b50905ec998
# ╠═10a37328-80bb-11eb-1899-c9abbcfaaaf0
# ╟─5e81dd2c-80b9-11eb-364e-63fa4ebd27da
# ╟─4f8dee16-80ba-11eb-1d0f-f3ecdc9ac671
# ╟─02299d2a-8010-11eb-024a-650734491c82
# ╟─90a78c38-80c2-11eb-3124-1dadb02c2806
# ╟─79a2aad2-7c2a-11eb-00b0-6bebd6658619
# ╟─968af7ba-7c2a-11eb-0aa4-6d51510020b1
# ╠═58e3898c-7c21-11eb-3364-9168fcedb22c
# ╟─cbf8bba4-7c26-11eb-3301-07233740f82a
# ╟─d1644026-7c27-11eb-17ae-67c414382cdc
# ╠═0ddf4fe6-7c28-11eb-2969-7b856ff408be
# ╟─32f4bc7c-80b9-11eb-0111-074cc68e892f
# ╠═728d5546-7a96-11eb-1757-93681c4e7955
# ╠═ecad4e78-7a96-11eb-2b54-e7f910e83835

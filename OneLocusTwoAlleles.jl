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

# ╔═╡ 8491344a-7c39-11eb-327e-df48b337d147
begin
	using Plots
	using PlutoUI
	using Measures
end

# ╔═╡ 0904e59a-7c3b-11eb-307d-6fd2049ea8bc
md"To use Eulers Method on an ODE with parameters you can either define global variables and use them within the function or to increase runtime pass all the parameters in a variable P to the vector field."

# ╔═╡ 8f2743fe-7c39-11eb-27d6-d527b21f46b8
begin
	function euler(f,t₀,tₑ,Δt,x₀,P)
		T = t₀:Δt:tₑ
		F = Vector{typeof(x₀)}(undef,length(T))

		F[1] = x₀
		
		for (n,tₙ) ∈ enumerate(@view T[2:end])
			F[n+1] = F[n] .+ f(tₙ,F[n],P) .* Δt
		end
		return F
	end
	
	#reshape result for plotting
	reshape_result(F) = [[f[i] for f ∈ F] for i ∈ 1:length(F[1])]
end;

# ╔═╡ 925558e8-7cc7-11eb-37a8-c98fb8a33355
md"""
---
### One Locus, Two Alleles

We describe a population where every individual is characterized by one locus only. There are two possible alleles 0 and 1, where 0 is the wildtype and 1 the mutated allele. Individuals who carry the mutated allel on both copies of the allele are thought to have a severe diseas, causing them to be excluded from the birth process. Individuals give birth and die at constant rates b, d. The carrying capacity Kₜ is time dependent and jumps at a certain point in time. Mutations appear at every birth at a constant rate μ from 0 to 1 only. The deterministic system can be described by the following system of non linear first order ODEs:

${\frac{dX_{00}}{dt} = b (1-\mu) \frac{(X_{00}+\frac{1}{2}X_{01})^2}{X_{00}+X_{01}} - (d + c\Sigma)X_{00} }$ 

${\frac{dX_{01}}{dt} = b \left((1-\mu) \frac{(X_{00}+\frac{1}{2}X_{01})X_{01}}{X_{00}+X_{01}} + \mu  \frac{(X_{00}+\frac{1}{2}X_{01})^2}{X_{00}+X_{01}} \right)- (d + c(\Sigma))X_{01} }$ 

${\frac{dX_{11}}{dt} = b \left((1-\mu) \frac{\frac{1}{4}X_{01}^2}{X_{00}+X_{01}} + \mu  \frac{(X_{00}+\frac{1}{2}X_{01})X_{01}}{X_{00}+X_{01}} + \mu^2  \frac{(X_{00}+\frac{1}{2}X_{01})^2}{X_{00}+X_{01}} \right)- (d + c(\Sigma))X_{01} }$ 

"""

# ╔═╡ 62311b66-7c3b-11eb-335c-c53b566fcb9f
begin	
	function OneLocus(t,X,P)
		X₀₀, X₀₁, X₁₁ = (x for x ∈ X)
		b, d, K, μ = (p for p ∈ P)
		
		N = sum(X)
		M = N - X₁₁
		D = d+((b-d)/K(t))*N
		
		return [
				(b * (1-μ)*(X₀₀+0.5*X₀₁)^2/M - D*X₀₀),
				(b * ((1-μ)*(X₀₀+0.5*X₀₁)*X₀₁/M + μ*(X₀₀+0.5*X₀₁)^2/M) - D*X₀₁),
				(b * ((1-μ)*0.25*(X₀₁)^2/M + μ*(X₀₀+0.5*X₀₁)*X₀₁/M + μ^2*(X₀₀+0.5*X₀₁)^2/M) - D*X₁₁)
				]
		end
end

# ╔═╡ 348cb7ba-7c40-11eb-1a80-4768b0c6b4d4
md"""
birth rate: $(@bind b Slider(0.0:0.1:1.0, show_value=true, default=1.0)) |
death rate: $(@bind d Slider(0.0:0.1:1.0, show_value=true, default=0.9)) |

initial population: $(@bind Kstart NumberField(1:10^7, default=500.0)) |
final population: $(@bind Kend NumberField(1:10^7, default=10_000.0)) |
burn in time: $(@bind tburn Slider(0:300,show_value=true,default=100))

mutation rate: $(@bind μ NumberField(0.0:0.0001:0.2, default=0.001))
"""

# ╔═╡ a19430e8-7c3d-11eb-0790-4b883ac9cd39
begin
	t_start = 0
	t_end = 500
	step = 0.01

	OL_0 = [Kstart,0.0,0.0]
	
	T = t_start:step:t_end
	K(t) = t <= tburn ? Kstart : Kend ;
	
	F = euler(OneLocus,t_start,t_end,step,OL_0,[b,d,K,μ])
	rF = reshape_result(F)
end;

# ╔═╡ d728c114-7c40-11eb-1fcb-a1ecbb31ae54
let
	p = plot(framestyle=:zerolines)
	
	plot!(T,rF,label=["(0,0)" "(0,1)" "(1,1)"])
	
	p
end	

# ╔═╡ 20a4f4dc-7c41-11eb-26a8-21e34fcfd8e4
begin
	mutation_load(X) = (X[2] + 2*X[3])/sum(X)
	ill_individual(X) = X[3]/sum(X)
end;

# ╔═╡ 619eab70-7c41-11eb-1fc0-fb3ec4e64220
let
	p_left = plot(framestyle=:zerolines,legend=:topleft,rightmargin=12mm,size=(680,400))
	
	plot!(p_left,T,mutation_load.(F),label="mutation load",color=:red)
	
	p_right = twinx()
	plot!(p_right,framestyle=:zerolines,legend=:topright,grid=false)
	plot!(p_right,T,ill_individual.(F),label="ill_individual",color=:orange)
	
	vline!(p_right,[T[end]],color=:black,label="")
	
end

# ╔═╡ Cell order:
# ╟─8491344a-7c39-11eb-327e-df48b337d147
# ╟─0904e59a-7c3b-11eb-307d-6fd2049ea8bc
# ╠═8f2743fe-7c39-11eb-27d6-d527b21f46b8
# ╟─925558e8-7cc7-11eb-37a8-c98fb8a33355
# ╠═62311b66-7c3b-11eb-335c-c53b566fcb9f
# ╟─348cb7ba-7c40-11eb-1a80-4768b0c6b4d4
# ╟─d728c114-7c40-11eb-1fcb-a1ecbb31ae54
# ╟─619eab70-7c41-11eb-1fc0-fb3ec4e64220
# ╠═a19430e8-7c3d-11eb-0790-4b883ac9cd39
# ╟─20a4f4dc-7c41-11eb-26a8-21e34fcfd8e4

### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 4b35680a-1310-11ed-29fa-3501f54a3932
begin
	using PlutoUI
	using Plots
	using Distributions
end

# ╔═╡ b970c566-bcb7-40bb-baed-598975411f14
md"Population size power $(@bind n Slider(1:10,default=5,show_value=true))"

# ╔═╡ 46f8a038-bb56-4d3d-b5fb-cfdf04afc067
md"Fitness coefficient s: $(@bind s Slider(0:0.01:1,default=0.02,show_value=true))"

# ╔═╡ 0acf9156-4c35-4b5c-954e-640849faecd4
md"Mutation rate per genome per generation $(@bind λ Slider(0.0:0.01:2.0,default=0.2,show_value=true))"

# ╔═╡ 5a198b4a-a549-470c-90f0-8203596e5ebb
md"Time: $(@bind Tend Slider(1000:1000:30_000,default=1000,show_value=true))"

# ╔═╡ 539de6af-7737-4e18-a307-7e62ac5dfb87
begin
	N = 10^n

	t_end = Tend+1
	slots = 50
end;

# ╔═╡ 632a74ee-71bc-43c7-a2e3-3350cc69a120
 md"""
 Inital distribution: $(@bind inital Select(["healthy" => "healthy", "equilibrium" => "equilibrium"]))
 """

# ╔═╡ bc230823-68f7-4633-9c90-d3dcc4de4877
md"#### Plots"

# ╔═╡ bbbfef59-0c77-414f-98f8-b3927d7a63bf
md"Time $(@bind t Slider(1:t_end,show_value=true,default=200)) Show Extinction Times? $(@bind showtimes CheckBox())"

# ╔═╡ e1ba4486-bd7d-431d-9aa3-ccd39dd7d932
md"show theoretical equilibrium? $(@bind showeq CheckBox())"

# ╔═╡ 046db000-4257-47bc-9e9c-d47c2bfd15bf
function mean_extinction_time(Fit,J,n_runs)
	#if for all runs this class did die out
	if Fit[end-1,J+1] == 0
		exs = [Fit[t,J+1]-Fit[t+1,J+1] for t in 1:size(Fit)[1]-1]
		return sum( (t-1) * p for (t,p) in enumerate(exs))/n_runs
	else 
		return Inf
	end
end

# ╔═╡ 501f9caf-f954-4678-9d26-84e177e97013
#Gives the relative frequency of populations where untile time t J fittest populations did die out for Pₙ(τᴶ>t) 
extinction_prob(Fit,J,t,n_runs) = (n_runs - Fit[t,J+1])/n_runs

# ╔═╡ f7774f00-e2a0-4014-a27a-5515b7b1a463
md"""
---
##### Model Functions
"""

# ╔═╡ 613db540-e773-413f-80ab-633113db3cff
begin
	p(k,X,s,Poi::Distribution,T) = sum(
		X[k+2-j]*(1-s)^(k+X[1]-j) * pdf(Poi,j-X[1]) / T
		for j in X[1]:k
			)
	p(k,X,s,Poi::Distribution) = p(k,X,s,Poi::Distribution,T(1,X,s))
	p(k,X,s,λ) = p(k,X,s,Poisson(λ))

	T(r,X,s) = sum(X[i]*(1-s)^((X[1]+i-2)*r) for i in 2:length(X))
end

# ╔═╡ 02a23239-61db-4e8a-965f-2388430d2ff5
function runHaigh(X₀,s,λ,t_end;ε=10^(-8))
	
	K = length(X₀)-1
	N = sum(view(X₀,2:K+1))
	Xs = Matrix{Integer}(undef,t_end,K+1)
	Xs[1,:] .= X₀
	ps = zeros(K)
	Poi = Poisson(λ)

	for t in 2:t_end
		T₁ = T(1,view(Xs,t-1,:),s)
		ps .= [p(i,view(Xs,t-1,:),s,Poi,T₁) for i in Xs[t-1,1]:Xs[t-1,1]+(K-1)]
		while ps[end] > ε
			Xs[t-1,1] += 1
			for i in 2:K
				Xs[t-1,i] = Xs[t-1,i+1]
			end
			popfirst!(ps)
			push!(ps,p(K+Xs[t-1,1]-1,view(Xs,t-1,:),s,Poi))
		end
		Xs[t,1] = Xs[t-1,1]
		Xs[t,2:K+1] .= rand(Multinomial(N,ps))
	end

	return Xs
end

# ╔═╡ d60f0c8e-2756-4ad7-a63e-0740c47999f1
function findextime(Xs,i)
	extime = findfirst(iszero,Xs[:,i])
	if extime == 1
		birthtime = findfirst(x->x>10,Xs[:,i])
		t_end = length(Xs[:,i])
		extime = findfirst(iszero,Xs[birthtime:t_end,i])
	end
	return extime
end

# ╔═╡ 82d91604-4a92-40eb-b293-2007aa284b23
allextimes(X,t,slot) = [findextime(X,i) for i in 2:findfirst(!iszero,X[t,2:slot+1])]

# ╔═╡ 141ac859-6600-467f-9be6-b2cdbf8f8b1e
begin
	eq(k,N,Poi::Distribution) = N*pdf(Poi,k)
	eq(k,N,θ) = eq(k,N,Poisson(θ))
end

# ╔═╡ 76726853-7b3c-475c-82c8-6844a06e1d81
begin
	load(X) = sum(number_mutations(X,i)*X[i] for i in 2:length(X))
	number_mutations(X,i) = X[1]+i-1
	mean_fitness(X,s) = sum(X[i]*(1-s)^number_mutations(X,i-1) for i in 2:length(X))/(sum(X)-X[1])
end

# ╔═╡ 0a49f851-4a58-4545-877a-f452bb2e47bf
function eq_initial(N,λ,s,slots)
	
	Poi = Poisson(λ/s)
	p = [pdf(Poi,k) for k in 0:slots-1]
	X₀ = vcat([0],round.(Integer,N*p))

	return X₀
end

# ╔═╡ 50165c8f-1fad-4e46-99f1-0ed242fd996f
begin
	if inital == "healthy"
		X₀ = vcat([0,N],zeros(Integer,slots-1))
	elseif inital == "equilibrium"
		X₀ = eq_initial(N,λ,s,slots)
	end
end;

# ╔═╡ 7010e148-d503-4ee9-96df-a0a520177eda
Xs = runHaigh(X₀,s,λ,t_end)

# ╔═╡ a973fe35-ffac-4d03-96a3-5e1f7db64add
ML = [load(Xs[i,:]) for i in 1:t_end];

# ╔═╡ 7c674ba9-1aa8-48d6-8201-7120dce75a92
extimes = allextimes(Xs,t_end,slots);

# ╔═╡ 7b53600c-db85-425a-81b6-cbd4e331637b
begin
	plot([mean_fitness(Xs[i,:],s) for i in 1:t_end],ylabel="Mean Fitness",xlabel="Time",legend=false)
	hline!([(1-s)^k * exp(-λ) for k in 0:length(extimes)])
end

# ╔═╡ 31952d34-df98-4369-bc1e-be11abd7c307
begin
	N_mean = 10^5
	s_mean = 0.02
	λ_mean = 0.2

	t_end_mean = 1_000+1
	slots_mean = 50
	n_runs = 50

	Xss = []
	Fit = zeros(Int,(t_end_mean,slots))
	
	for run in 1:n_runs
		Xs = runHaigh(eq_initial(N_mean,λ_mean,s_mean,slots_mean),s_mean,λ_mean,t_end_mean)
		
		extimes = allextimes(Xs,t_end_mean,slots_mean)
		extimes = vcat([0],extimes,[t_end_mean-1])
		
		for i in 1:length(extimes)-1
			Fit[extimes[i]+1:extimes[i+1],i] .+= 1
		end

		push!(Xss,Xs)
	end
	
	Xs_mean = mean(Xss)
end;

# ╔═╡ 69a257d6-7689-46a6-b132-354b29290c77
md"##### Mean over $(n_runs) Simulations"

# ╔═╡ c485f4a0-ce6e-4eb8-b968-b737bd347bea
mean_extinction_time(Fit,0,n_runs)

# ╔═╡ bb65123e-7815-4a5a-923b-567e59af475a
function fittestclass(extimes,t) 
	exs = findfirst([ex > t for ex in extimes])
	if isnothing(exs)
		return length(extimes)
	else
		return exs-1
	end
end

# ╔═╡ 5d0ee2af-d22b-49a8-b071-2e631c181e00
begin
	pt = plot(ML,xlabel="Mutation Load",legend=false)
	vline!(pt,[t])
	if showtimes && !isempty(extimes)
		vline!(pt,extimes)
	end
	br = bar(Xs[t,2:slots+1],ylim=(0,1.5*10^(n-1)),xlabel="Classes",legend=false)
	fittest = fittestclass(extimes,t)
	bar!(br,[eq(k-fittest,N,Poisson(λ/s)) for k in 0:slots],color=:orange,alpha=0.5)
	plot(pt,br,layout=(2,1),size=(900,1200))
end

# ╔═╡ cee0add3-1172-4247-a2d0-c702e87df32a
function calculate_eqs(Fit,t,N,λ,s,slots,n_runs)
	EQ = zeros(slots)
	P=Poisson(λ/s)
	for (J,p) in enumerate(Fit[t,:])
		if p > 0 
			EQ .+= p/n_runs .* [eq(k-(J-1),N,P) for k in 0:slots-1]
		end
	end
	return EQ
end

# ╔═╡ 6a093cb7-6834-4ef6-86cb-7abc92327ff3
begin
	Ts = [0,200,400,600,800,1000]
	pl = plot()

	for t in Ts
		plot!(pl, 0:slots_mean-1,Xs_mean[t+1,2:slots_mean+1], 
			label="t=$t",
			xlim=(0,30),
			yticks=([0,5*10^3,10*10^3,15*10^3],[0,5,10,15]),
			ylim=(0,15*10^3),
			ylabel = "Frequency x 10⁻³",
			xlabel = "Classes",
			palette=palette(:Blues,length(Ts))
		)
	end

	if showeq
		plot!(pl,
			[calculate_eqs(Fit,t+1,N_mean,λ_mean,s_mean,slots_mean,n_runs) for t in Ts],
			palette=palette(:Reds,length(Ts)),
			alpha = 0.75,
			label = ""
		)
	end
	
	pl
	#savefig("HaighPlot.pdf")
end

# ╔═╡ 50acd1e1-d03f-4cbb-a1d3-7d8c55fdc0b6
begin
	plot(0:slots_mean-1,
		[calculate_eqs(Fit,t+1,N_mean,λ_mean,s_mean,slots_mean,n_runs) for t in Ts],
		label=reshape(["t=$t" for t in Ts],(1,length(Ts))),
		xlim=(0,30),
		yticks=([0,5*10^3,10*10^3,15*10^3],[0,5,10,15]),
		ylim=(0,15*10^3),
		ylabel = "Frequency x 10⁻³",
		xlabel = "Classes",
		palette=palette(:Blues,length(Ts)),
		title = "Theoretical Equilibrium Sizes"
	)
end

# ╔═╡ 23f0b75e-6f0d-47de-93f4-b7fd992c1211
function recurrence(f,init,n_end)
	res = [f(init)]
	for i ∈ 1:n_end-1
		push!(res,f(res[end]))
	end
	return res
end

# ╔═╡ e25b2e0a-94cb-49c2-a069-503f0407f68e
#some important functions
begin
	#probability generating function of Poi(m)
	f(x,m) = exp(m*(x-1))
	#mean of offspring distribution for poisson branching process
	m(s,J) = 1/(1-s)^J
	#equilibrium population size of fittest class
	n₀(N,λ,s) = N*exp(-λ/s)
	#deterministic time to reestablish equilibrium
	d(λ,s) = -log(λ/s)/log(1-s)
end

# ╔═╡ 53f4aff2-cfda-4408-b57f-40805a82d25d
a = recurrence(x->f(x,m(s_mean,0)),0,1000) .^ n₀(N_mean,λ_mean,s_mean)

# ╔═╡ 3d0fc2ac-5d08-4592-872b-36ee2d3a01d1
plot([[extinction_prob(Fit,0,t,n_runs) for t in 1:1000],a])

# ╔═╡ 3a514225-a4a5-4e94-a376-b002151a536c
n₀(N_mean,λ_mean,s_mean)

# ╔═╡ 645c06f2-67c2-49fa-b4b0-551a0d7c98bd
d(λ_mean,s_mean)

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Distributions = "31c24e10-a181-5473-b8eb-7969acd0382f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Distributions = "~0.25.71"
Plots = "~1.24.3"
PlutoUI = "~0.7.40"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

[[AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[ArgTools]]
git-tree-sha1 = "bdf73eec6a88885256f282d48eafcad25d7de494"
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[Artifacts]]
deps = ["Pkg"]
git-tree-sha1 = "c30985d8821e0cd73870b17b0ed0ce6dc44cb744"
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"
version = "1.3.0"

[[Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c3598e525718abcc440f69cc6d5f60dda0a1b61e"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.6+5"

[[Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "e2f47f6d8337369411569fd45ae5753ca10394c6"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.0+6"

[[Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "32ad4ece064a61855a35bdc34e3da0b495e01169"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.12.2"

[[ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[Compat]]
deps = ["Base64", "Dates", "DelimitedFiles", "Distributed", "InteractiveUtils", "LibGit2", "Libdl", "LinearAlgebra", "Markdown", "Mmap", "Pkg", "Printf", "REPL", "Random", "SHA", "Serialization", "SharedArrays", "Sockets", "SparseArrays", "Statistics", "Test", "UUIDs", "Unicode"]
git-tree-sha1 = "78bee250c6826e1cf805a88b7f1e86025275d208"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "3.46.0"

[[CompilerSupportLibraries_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "8e695f735fca77e9708e795eda62afdb869cbb70"
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.3.4+0"

[[Contour]]
deps = ["StaticArrays"]
git-tree-sha1 = "9f02045d934dc030edad45944ea80dbd1f0ebea7"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.5.7"

[[DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[DensityInterface]]
deps = ["InverseFunctions", "Test"]
git-tree-sha1 = "80c3e8639e3353e5d2912fb3a1916b8455e2494b"
uuid = "b429d917-457f-4dbc-8f4c-0cc954292b1d"
version = "0.4.0"

[[Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[Distributions]]
deps = ["ChainRulesCore", "DensityInterface", "FillArrays", "LinearAlgebra", "PDMats", "Printf", "QuadGK", "Random", "SparseArrays", "SpecialFunctions", "Statistics", "StatsBase", "StatsFuns", "Test"]
git-tree-sha1 = "ee407ce31ab2f1bacadc3bd987e96de17e00aed3"
uuid = "31c24e10-a181-5473-b8eb-7969acd0382f"
version = "0.25.71"

[[DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
git-tree-sha1 = "39e99578597b4b1660b63cdabd5224ba53e3e71a"
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[DualNumbers]]
deps = ["Calculus", "NaNMath", "SpecialFunctions"]
git-tree-sha1 = "5837a837389fccf076445fce071c8ddaea35a566"
uuid = "fa6b7ba4-c1ee-5f82-b5fc-ecf0adba8f74"
version = "0.6.8"

[[EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "92d8f9f208637e8d2d28c664051a00569c01493d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.1.5+1"

[[Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1402e52fcda25064f51c77a9655ce8680b76acf0"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.2.7+6"

[[FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "LibVPX_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "3cc57ad0a213808473eafef4845a74766242e05f"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.3.1+4"

[[FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[FillArrays]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "Statistics"]
git-tree-sha1 = "deed294cde3de20ae0b2e0355a6c4e1c6a5ceffc"
uuid = "1a297f60-69ca-5386-bcde-b61e274b549b"
version = "0.12.8"

[[FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "35895cf184ceaab11fd778b4590144034a167a2f"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.1+14"

[[Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "cbd58c9deb1d304f5a245a0b7eb841a2560cfec6"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.1+5"

[[FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0d20aed5b14dd4c9a2453c1b601d08e1149679cc"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.5+6"

[[GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "a199aefead29c3c2638c3571a9993b564109d45a"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.4+0"

[[GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "30f2b340c2fff8410d89bfcdc9c0a6dd661ac5f7"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.62.1"

[[GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "d59e8320c2747553788e4fc42231489cc602fa50"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.58.1+0"

[[GeometryBasics]]
deps = ["EarCut_jll", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "83ea630384a13fc4f002b77690bc0afeb4255ac9"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.2"

[[Gettext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "8c14294a079216000a0bdca5ec5a447f073ddc9d"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.20.1+7"

[[Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "04690cc5008b38ecbdfede949220bc7d9ba26397"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.59.0+4"

[[Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[HTTP]]
deps = ["Base64", "Dates", "IniFile", "Logging", "MbedTLS", "NetworkOptions", "Sockets", "URIs"]
git-tree-sha1 = "0fa77022fe4b511826b39c894c90daf5fce3334a"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "0.9.17"

[[HypergeometricFunctions]]
deps = ["DualNumbers", "LinearAlgebra", "OpenLibm_jll", "SpecialFunctions", "Test"]
git-tree-sha1 = "709d864e3ed6e3545230601f94e11ebc65994641"
uuid = "34004b35-14d8-5ef3-9330-4cdb6864b03a"
version = "0.3.11"

[[Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9aff0587d9603ea0de2c6f6300d9f9492bbefbd3"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.0.1+3"

[[LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "df381151e871f41ee86cee4f5f6fd598b8a68826"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.0+3"

[[LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f128cd6cd05ffd6d3df0523ed99b90ff6f9b349a"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.0+3"

[[LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
git-tree-sha1 = "cdbe7465ab7b52358804713a53c7fe1dac3f8a3f"
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[LibCURL_jll]]
deps = ["LibSSH2_jll", "Libdl", "MbedTLS_jll", "Pkg", "Zlib_jll", "nghttp2_jll"]
git-tree-sha1 = "897d962c20031e6012bba7b3dcb7a667170dad17"
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.70.0+2"

[[LibGit2]]
deps = ["Printf"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[LibSSH2_jll]]
deps = ["Libdl", "MbedTLS_jll", "Pkg"]
git-tree-sha1 = "717705533148132e5466f2924b9a3657b16158e8"
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.9.0+3"

[[LibVPX_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "85fcc80c3052be96619affa2fe2e6d2da3908e11"
uuid = "dd192d2f-8180-539f-9fb4-cc70b1dcf69a"
version = "1.9.0+1"

[[Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a2cd088a88c0d37eef7d209fd3d8712febce0d90"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.1+4"

[[Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "b391a18ab1170a2e568f9fb8d83bc7c780cb9999"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.5+4"

[[Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ec7f2e8ad5c9fa99fc773376cdbc86d9a5a23cb7"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.36.0+3"

[[Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "cba7b560fcc00f8cd770fa85a498cbc1d63ff618"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.0+8"

[[Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51ad0c01c94c1ce48d5cad629425035ad030bfd5"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.34.0+3"

[[Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "291dd857901f94d683973cdf679984cdf73b56d0"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.1.0+2"

[[Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f879ae9edbaa2c74c922e8b85bb83cc84ea1450b"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.34.0+7"

[[LinearAlgebra]]
deps = ["Libdl"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "Random", "Sockets"]
git-tree-sha1 = "1c38e51c3d08ef2278062ebceade0e46cefc96fe"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.0.3"

[[MbedTLS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0eef589dd1c26a3ac9d753fe1a8bcad63f956fa6"
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.16.8+1"

[[Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[MozillaCACerts_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e49d1b24d15a6995a5781b6670dc4a9d5c5aef8d"
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.3.29+0"

[[NaNMath]]
git-tree-sha1 = "b086b7ea07f8e38cf122f5016af580881ac914fe"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "0.3.7"

[[NetworkOptions]]
git-tree-sha1 = "ed3157f48a05543cce9b241e1f2815f7e843d96e"
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "a42c0f138b9ebe8b58eba2271c5053773bde52d0"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.4+2"

[[OpenLibm_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "d22054f66695fe580009c09e765175cbf7f13031"
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.7.1+0"

[[OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "71bbbc616a1d710879f5a1021bcba65ffba6ce58"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.1+6"

[[OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9db77584158d0ab52307f8c04f8e7c08ca76b5b3"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.3+4"

[[Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f9d57f4126c39565e05a2b0264df99f497fc6f37"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.1+3"

[[OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b556ad51dceefdbf30e86ffa8f528b73c7df2bb"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.42.0+4"

[[PDMats]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "cf494dca75a69712a72b80bc48f59dcf3dea63ec"
uuid = "90014a1f-27ba-587c-ab20-58faa44d9150"
version = "0.11.16"

[[Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6a20a83c1ae86416f0a5de605eaea08a552844a3"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.0+0"

[[Pkg]]
deps = ["Dates", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "UUIDs"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"

[[PlotThemes]]
deps = ["PlotUtils", "Requires", "Statistics"]
git-tree-sha1 = "a3a964ce9dc7898193536002a6dd892b1b5a6f1d"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "2.0.1"

[[PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "b084324b4af5a438cd63619fd006614b3b20b87b"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.0.15"

[[Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun"]
git-tree-sha1 = "d73736030a094e8d24fdf3629ae980217bf1d59d"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.24.3"

[[PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "a602d7b0babfca89005da04d89223b867b55319f"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.40"

[[Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "16626cfabbf7206d60d84f2bf4725af7b37d4a77"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.2+0"

[[QuadGK]]
deps = ["DataStructures", "LinearAlgebra"]
git-tree-sha1 = "3c009334f45dfd546a16a57960a821a1a023d241"
uuid = "1fd47b50-473d-5c70-9696-f719f8f3bcdc"
version = "2.5.0"

[[REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "7ad0dfa8d03b7bcf8c597f59f5292801730c55b8"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.4.1"

[[Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[Rmath]]
deps = ["Random", "Rmath_jll"]
git-tree-sha1 = "86c5647b565873641538d8f812c04e4c9dbeb370"
uuid = "79098fc4-a85e-5d69-aa6a-4863f24498fa"
version = "0.6.1"

[[Rmath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "1b7bf41258f6c5c9c31df8c1ba34c1fc88674957"
uuid = "f50d1b31-88e8-58de-be2c-1cc44531875f"
version = "0.2.2+2"

[[SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[StaticArrays]]
deps = ["LinearAlgebra", "Random", "Statistics"]
git-tree-sha1 = "74eaf352c0cef1e32ce7394bcc359d9199a28cf7"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.3.6"

[[Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[StatsFuns]]
deps = ["ChainRulesCore", "HypergeometricFunctions", "InverseFunctions", "IrrationalConstants", "LogExpFunctions", "Reexport", "Rmath", "SpecialFunctions"]
git-tree-sha1 = "5783b877201a82fc0014cbf381e7e6eb130473a4"
uuid = "4c63d2b9-4356-54db-8cca-17b64c39e42c"
version = "1.0.1"

[[StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArrays", "Tables"]
git-tree-sha1 = "ec47fb6069c57f1cee2f67541bf8f23415146de7"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.11"

[[SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[TOML]]
deps = ["Dates"]
git-tree-sha1 = "44aaac2d2aec4a850302f9aa69127c74f0c3787e"
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[Test]]
deps = ["Distributed", "InteractiveUtils", "Logging", "Random"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[URIs]]
git-tree-sha1 = "97bbe755a53fe859669cd907f2d96aee8d2c1355"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.3.0"

[[UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "dc643a9b774da1c2781413fd7b6dcd2c56bb8056"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.17.0+4"

[[Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "be0db24f70aae7e2b89f2f3092e93b8606d659a6"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.10+3"

[[XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "2b3eac39df218762d2d005702d601cd44c997497"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.33+4"

[[Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[Zlib_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "320228915c8debb12cb434c59057290f0834dbf6"
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.11+18"

[[Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "2c1332c54931e83f8f94d310fa447fd743e8d600"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.4.8+0"

[[libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "acc685bcf777b2202a904cdcb49ad34c2fa1880c"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.14.0+4"

[[libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7a5780a0d9c6864184b3a2eeeb833a0c871f00ab"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "0.1.6+4"

[[libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "6abbc424248097d69c0c87ba50fcb0753f93e0ee"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.37+6"

[[libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "fa14ac25af7a4b8a7f61b287a124df7aab601bcd"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.6+6"

[[nghttp2_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "8e2c44ab4d49ad9518f359ed8b62f83ba8beede4"
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.40.0+2"

[[x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "d713c1ce4deac133e3334ee12f4adff07f81778f"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2020.7.14+2"

[[x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "487da2f8f2f0c8ee0e83f39d13037d6bbf0a45ab"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.0.0+3"

[[xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "ece2350174195bb31de1a63bea3a41ae1aa593b6"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "0.9.1+5"
"""

# ╔═╡ Cell order:
# ╠═50165c8f-1fad-4e46-99f1-0ed242fd996f
# ╠═7010e148-d503-4ee9-96df-a0a520177eda
# ╠═a973fe35-ffac-4d03-96a3-5e1f7db64add
# ╠═7c674ba9-1aa8-48d6-8201-7120dce75a92
# ╠═539de6af-7737-4e18-a307-7e62ac5dfb87
# ╟─b970c566-bcb7-40bb-baed-598975411f14
# ╟─46f8a038-bb56-4d3d-b5fb-cfdf04afc067
# ╟─0acf9156-4c35-4b5c-954e-640849faecd4
# ╟─5a198b4a-a549-470c-90f0-8203596e5ebb
# ╟─632a74ee-71bc-43c7-a2e3-3350cc69a120
# ╟─bc230823-68f7-4633-9c90-d3dcc4de4877
# ╟─bbbfef59-0c77-414f-98f8-b3927d7a63bf
# ╠═5d0ee2af-d22b-49a8-b071-2e631c181e00
# ╟─7b53600c-db85-425a-81b6-cbd4e331637b
# ╟─69a257d6-7689-46a6-b132-354b29290c77
# ╟─e1ba4486-bd7d-431d-9aa3-ccd39dd7d932
# ╟─6a093cb7-6834-4ef6-86cb-7abc92327ff3
# ╟─50acd1e1-d03f-4cbb-a1d3-7d8c55fdc0b6
# ╠═31952d34-df98-4369-bc1e-be11abd7c307
# ╠═3d0fc2ac-5d08-4592-872b-36ee2d3a01d1
# ╠═53f4aff2-cfda-4408-b57f-40805a82d25d
# ╠═c485f4a0-ce6e-4eb8-b968-b737bd347bea
# ╠═3a514225-a4a5-4e94-a376-b002151a536c
# ╠═645c06f2-67c2-49fa-b4b0-551a0d7c98bd
# ╟─f7774f00-e2a0-4014-a27a-5515b7b1a463
# ╠═4b35680a-1310-11ed-29fa-3501f54a3932
# ╠═501f9caf-f954-4678-9d26-84e177e97013
# ╠═046db000-4257-47bc-9e9c-d47c2bfd15bf
# ╠═613db540-e773-413f-80ab-633113db3cff
# ╠═02a23239-61db-4e8a-965f-2388430d2ff5
# ╠═d60f0c8e-2756-4ad7-a63e-0740c47999f1
# ╠═82d91604-4a92-40eb-b293-2007aa284b23
# ╠═141ac859-6600-467f-9be6-b2cdbf8f8b1e
# ╠═76726853-7b3c-475c-82c8-6844a06e1d81
# ╠═0a49f851-4a58-4545-877a-f452bb2e47bf
# ╠═bb65123e-7815-4a5a-923b-567e59af475a
# ╠═cee0add3-1172-4247-a2d0-c702e87df32a
# ╠═23f0b75e-6f0d-47de-93f4-b7fd992c1211
# ╠═e25b2e0a-94cb-49c2-a069-503f0407f68e
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002

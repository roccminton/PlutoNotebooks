### A Pluto.jl notebook ###
# v0.12.12

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

# ╔═╡ 38ad5806-f5c2-11eb-1a4b-757cb20b199c
begin
	using Plots
	using PlutoUI
	using LinearAlgebra
	using DifferentialEquations
	using Statistics
	using CSV
	using DataFrames
end

# ╔═╡ a9dcf66c-00e1-11ec-2489-85a8bf50529f
Ns = 1:4;

# ╔═╡ e002f598-fcd6-11eb-0693-cd37b6731ddb
md"N $(@bind N Slider(Ns, show_value = true, default=2)) "

# ╔═╡ def9fb92-00e6-11ec-0ce0-8bf38dfba80b
dnis = range(0.0,stop=2*N,length=100);

# ╔═╡ e2864ca2-00eb-11ec-15af-13e2c905c2d9
md"De novo influx $(@bind dni Slider(dnis, show_value=true, default=0.0))"

# ╔═╡ edf4f43a-01c7-11ec-3f89-87267d95ba52
md"---
Produce Animations for R? $(@bind anims CheckBox(default=false))"

# ╔═╡ 634cae54-01c7-11ec-34df-9df3b4d9558d
md"--- 
Plot Equilibrium for all N? $(@bind plot_eq CheckBox(default=false))"

# ╔═╡ ff76b4b2-04e4-11ec-1e43-59b314d25a53
md"---"

# ╔═╡ c6b86160-04e5-11ec-050c-711567fe7ca8
csvN = [
	CSV.File("DataInput/1Gene_dni_eq_mutationload_prevalence.csv"),
	CSV.File("DataInput/2Genes_dni_eq_mutationload_prevalence.csv")
	]
	

# ╔═╡ 5ca4273e-01c8-11ec-078a-5798323bc39d
md"---"

# ╔═╡ b1a05bb4-00f0-11ec-1a54-13b82c896a6c
begin
	function plot_dni_bar(dni,dnis)
			p =plot(
				yaxis=false,
				legend=false,
				ylim=(0,2),
				ticks=false,
				grid=false,
				showaxis=false,
				xlabel="De novo influx"
			)
	
			plot!(p,dnis,fill(1.0,length(dnis)),color=:black)
			scatter!(p,[0.0,dnis[end]],[1.0,1.0],m = (1, :vline, 10),color=:black)
			scatter!(p,[dni],[1.35],m = (1, :vline, 7),color=:black)
	
			scatter!(p,[dni],[1.0],m = (1, :red, 12))
	
			annotate!(p,[
					(0,0.2,"0"),
					(dnis[end],0.2,"$(ceil(Int,dnis[end]))"),
					(dni,1.8,"$(round(dni,digits=2))")
					])
		
			return p
	end
	
	function plot_dni_bar(i)
			p =plot(
				yaxis=false,
				legend=false,
				ylim=(0,2),
				ticks=false,
				grid=false,
				showaxis=false,
				xlabel="De novo influx"
			)
			#black line
			plot!(p,range(0.0,stop=1.0,length=100),fill(1.0,100),color=:black)
			#start and end bar
			scatter!(p,[0.0,1.0],[1.0,1.0],m = (1, :vline, 10),color=:black)
			#red dot
			scatter!(p,[i/100],[1.0],m = (1, :red, 12))
	
			annotate!(p,[
					(0.0,0.2,"0"),
					(1.0,0.2,"2N")
					])
		
			return p
	end
end

# ╔═╡ 41c242fe-ff58-11eb-01b5-81fac3c399cd
begin
	L_guess(dni,N) = (2^N*dni)^(1/2)
	I_guess(dni,N) = dni/(2^N)
end

# ╔═╡ 18afc95a-fcca-11eb-2c88-792999cb5897
function R(μ)
	ν = 1-μ
	return [
		ν^2 ν^2/2 (ν^2)/4
		2*μ*ν (μ*ν+ν/2) (μ*ν/2+ν/2)
		μ^2 ((μ^2)/2+μ/2) ((μ^2)/4+μ/2+1/4)
	]
end

# ╔═╡ a445432e-fce7-11eb-33f3-b71379b6eddb
begin
	#Number of possible matings
	number_of_columns(N) = 2^(N-1) * (2^N + 1)
	#Numbere of all Configurations
	number_of_rows(N) = 3^N
end

# ╔═╡ c4fdb038-00f7-11ec-3011-c9183ec5fe36
function plot_histogram(r,N,fbins,ymax,xmax)
		h = histogram(
			vcat([r[i,:] for i ∈ 1:size(r)[1]]...),
			bins=range(
				0.0,
				stop=1,
				length=ceil(Integer,number_of_columns(N)*number_of_rows(N)*fbins)
				),
			ylim=(0,ymax),
			xlim=(0,xmax),
			ylabel = "Number of Outcomes",
			xlabel = "Probability of particular outcome",
			legend = false,
			title="N=$N"
			)
end

# ╔═╡ 73f9e942-01cc-11ec-3482-d58071c89d6e
begin
	LoadN(x,w) = sum(x.*w)
	LoadN_healthy(x,w,h) = sum(x.*w.*h)
	LoadN_ill(x,w,h) = sum(x.*w.*abs.(h.-1))
	PrevN(x,h) = 1-sum(x.*h)
end

# ╔═╡ 5e6f439c-fe8e-11eb-1e33-2f7be0eaadff
begin
	ham_weight(bitstring) = sum( c == '0' for c ∈ bitstring)
	Zeros(N) = ham_weight.([string(i,base=2,pad=N) for i ∈ 0:2^N-1])
end

# ╔═╡ 71df4b56-fe73-11eb-3559-55a35d29f37c
begin
	Prev(x) = 1-sum(x)
	
	Load(x,N,z) = 2*N - sum(z[i]*x[i] for i ∈ 1:2^N)
	
	function Load(x)
		N = Int(log(2,size(x,1)))
		z = Zeros(N)
		return Load(x,N,z)
	end
end

# ╔═╡ 95801cc2-ff37-11eb-28c2-b16f7958283a
function limit(arr,acc=10^(-3))
	last = length(arr)
	i = findlast([abs(arr[i-1]-arr[i]) > acc for i in 2:last])
	isnothing(i) ? (return arr[end]) : (return mean(view(arr,i:last)))
end

# ╔═╡ c0d9a06c-fe78-11eb-3776-1723db5d3551
mutation_rate(dni,N)=dni/(2*N)

# ╔═╡ 6b6402f6-fcde-11eb-0a18-63cc8c83fae8
begin
	function convert_index_to_number(i,j,m)
		# m = 2^N the Number healty configurations
		return Int((m/2)*(m+1)-1/2*(m-i)*(m-i+1)-(m-j))
	end
	convert_index_to_number((i,j),m) = convert_index_to_number(i,j,m)
end

# ╔═╡ 8e4d6f14-fdef-11eb-1e04-278b64849bff
function diff_eq(x,hR,t)
	m = size(x,1) #equals 2^N
	Σ = sum(x)^2
	return [
			sum(
				sum(
					(x[i]*x[j]/Σ)*hR[k,convert_index_to_number(i,j,m)]
				for j ∈ 1:m) 
			for i ∈ 1:m)
		for k in 1:m] .- x
end

# ╔═╡ 670bcc9c-01d0-11ec-16a2-417c2d8d1924
function diff_eq_all(x,RHImn,t)
	r,H,I,m,n = RHImn
	Z = 1/(sum(x .* H)^2)
	return [
			Z * sum( H[i] &&
				sum( H[j] &&
					x[i]*x[j]*r[k,convert_index_to_number(I[i],I[j],n)]
				for j ∈ i:m) 
			for i ∈ 1:m)
		for k in 1:m] .- x
end

# ╔═╡ 27cbec3a-01d3-11ec-01c1-a7495311e4af
index_in_healthy_only(i,H) = sum(view(H,1:i))

# ╔═╡ f7c16cac-fcde-11eb-1a1e-37cbc0e31289
function convert_number_to_index(n,N)
	#Number of healthy configurations
	M = 2^N
	#Number of possible matings
	N = number_of_rows(N)
	acc = 0
	for (i,m) ∈ enumerate(M:-1:1)
		n <= m + acc ? (return Int(i),Int(n-acc+i-1)) : acc += m
	end
	return Int(N), Int(N)
end

# ╔═╡ abb776ea-00cf-11ec-3d9a-dbb357f9883e
begin
	bit_string(i,N) = string(i,base=2,pad=N)
	ter_string(i,N) = string(i,base=3,pad=N)
	convert_index_to_bin_string(i,j,N) = "[" * bit_string(i-1,N) * "]x[" *  bit_string(j-1,N) * "]"
	convert_index_to_bin_string((i,j),N) = convert_index_to_bin_string(i,j,N)
	
	xticks_heatmap(N) = (1:number_of_columns(N),convert_index_to_bin_string.(
		[convert_number_to_index(i,N) for i in 1:number_of_columns(N)]
		,N))
	yticks_heatmap(N) = (1:number_of_rows(N),reverse(["[" * ter_string(i-1,N) * "]" for i in 1:number_of_rows(N)]))
end

# ╔═╡ 0eb45f48-ff58-11eb-2704-4f5dad08241b
begin
	plot_heatmap(r,N,ticks::Bool,dni,cmax) = heatmap(
		reverse(r,dims=1), 
		color = :thermal, 
		xticks = (ticks ? xticks_heatmap(N) : false),
		yticks = (ticks ? yticks_heatmap(N) : false),
		xmirror = true,
		clims = (0.0,cmax),
		title="N=$N"
		) 
	
	function plot_heatmap(r,N,dni,dnis,cmax=1.0)
		plot(
			(N < 3 ? plot_heatmap(r,N,true,dni,cmax) : plot_heatmap(r,N,false,dni,cmax)),
			plot_dni_bar(dni,dnis),
			layout = grid(2, 1, heights=[0.79,0.2]),
			size=(700,500)
			)
	end
end

# ╔═╡ bedb15aa-fced-11eb-0f8c-6b0f1063b26f
begin
	get_binary_digit(num,digit) = floor(Int,mod(num/(2^(digit-1)),2))
	get_ternary_digit(num,digit) = floor(Int,mod(num/(3^(digit-1)),3))
	get_Nary_digit(num,digit,N) = floor(Int,mod(num/(N^(digit-1)),N))
end

# ╔═╡ 685cf83e-fce7-11eb-38ee-cb456fdbc26d
function generate_R(r,N)
	n, m = number_of_rows(N), number_of_columns(N)
	R = Matrix(undef, n, m)
	for col in 1:m
		i,j = convert_number_to_index(col,N)
		binsum = [get_binary_digit(i-1,k) + get_binary_digit(j-1,k) + 1 for k in 1:N]
		for row in 1:n
			tercol = [get_ternary_digit(row-1,k) + 1 for k in 1:N]
			R[row,col] = prod(r[a,b] for (a,b) in zip(tercol,binsum))
		end
	end
	return R
end

# ╔═╡ f3d2bb1c-00d0-11ec-3fc6-b7db47f9232e
R_all(N) = [generate_R(R(mutation_rate(dni,N)),N) for dni in range(0.0,stop=2*N,length=100)]

# ╔═╡ b1cb8c58-00e1-11ec-060f-7bc5f20723c6
All_R = [R_all(n) for n ∈ Ns];

# ╔═╡ 45ee92e0-00ee-11ec-2060-c791823a5884
let
	if anims
		anim = @animate for dni ∈ dnis

			r = All_R[N][round(Integer,dni*(length(dnis)-1)/(2*N) + 1)]

			plot(
				plot_histogram(r,N,0.1,1000,0.5),
				plot_dni_bar(dni,dnis),
				layout=grid(2, 1, heights=[0.79,0.2]),
				size=(700,500)
			)

		end
		gif(anim, "PlotOutput/MatingGridHist_N$N.gif", fps = 15)
	end
end

# ╔═╡ 478f641c-00ec-11ec-1f4a-cd703ab756d6
let
	if anims
		anim = @animate for dni ∈ dnis
			plot_heatmap(
				All_R[N][round(Integer,dni*(length(dnis)-1)/(2*N) + 1)],
				N,
				dni,
				dnis,
				2.0^(min(2-N,0))
			)
		end
		gif(anim, "PlotOutput/MatingGrid_N$N.gif", fps = 15)
	end
end

# ╔═╡ df0f60ec-fdf1-11eb-23e0-9f2747d9d223
begin
	ternary_list(N) = reverse.(Iterators.product(fill(0:2,N)...))[:]
	not_ill(tup) = 2 ∉ tup
	healthy_config_positions(N) = [not_ill(num) for num in ternary_list(N)]
end

# ╔═╡ acdef9b4-fdf3-11eb-0aaa-6de346bde1cd
function only_healthy(R,N)
	P = Matrix{Float64}(undef,2^N,number_of_columns(N))
	H = healthy_config_positions(N)
	j = 1
	for i ∈ 1:number_of_rows(N)
		if H[i] 
			P[j,:] = R[i,:]
			j += 1
		end
	end
	return P
end

# ╔═╡ 7d72ff6e-03f4-11ec-1b98-8f66c650f673
parameter(r,N) = [
	r,
	healthy_config_positions(N),
	cumsum(healthy_config_positions(N)),
	3^N,
	2^N
]

# ╔═╡ 6811d1c2-fe6a-11eb-05f8-151a60d2bd0a
begin
	r = All_R[N][round(Integer,dni*(length(dnis)-1)/(2*N) + 1)]
	
	x₀ = zeros(3^N)
	x₀[1] = 0.99
	x₀[2] = 0.01
	tspan = (0.0,25.0)
	H = healthy_config_positions(N)
	
	prob = ODEProblem(diff_eq_all, x₀, tspan, parameter(r,N))
	sol = solve(prob)
end;

# ╔═╡ e2278904-37e3-11ec-07e6-578d2e0e155c
begin
	plot(sol,label=["X⁰" "X¹" "X²"])
	savefig("PlotOutput/OneSimulation_N=1_mu=025.pdf")
end

# ╔═╡ 19aeceb6-0bcd-11ec-2237-a19898d1df26
plot([[u[k] for u in sol.u] for k in 1:length(sol.u[1])])

# ╔═╡ 2dd78c3a-00e6-11ec-33a2-07189bcd1803
plot_heatmap(r,N,dni,dnis)

# ╔═╡ a01907ac-03ee-11ec-3f15-81f0a341a1bb
begin
	char_to_int(c) = Int(c)-48
	string_to_int(s) = sum(char_to_int(c) for c ∈ s)
	Load_weights(N) = string_to_int.([string(i-1,base=3,pad=N) for i=1:3^N])
end

# ╔═╡ 91604328-ff38-11eb-0441-595f7d6f456f
begin
	if plot_eq
		eq_Load = Matrix{Float64}(undef,length(Ns),100)
		eq_Load_ill = Matrix{Float64}(undef,length(Ns),100)
		eq_Prev = Matrix{Float64}(undef,length(Ns),100)
		Sol = Matrix{Any}(undef,length(Ns),100)

		for (col,N) ∈ enumerate(Ns)

			dnis = range(0.0,stop=2*N,length=100)
			H = healthy_config_positions(N)
			W = Load_weights(N)

			for (row,dni) in enumerate(dnis)

				r = All_R[N][round(Integer,(dni*99)/(2*N) + 1)]

				x₀ = zeros(3^N)
				x₀[1] = 0.99
				x₀[2] = 0.01
				tspan = (0.0,100.0)

				prob = ODEProblem(diff_eq_all, x₀, tspan, parameter(r,N))
				sol = solve(prob)
				
				Sol[col,row] = sol

				eq_Load[col,row] = limit(
					map(x->LoadN(x,W),sol.u)
				)
				eq_Load_ill[col,row] = limit(
					map(x->LoadN_ill(x,W,H),sol.u)
				)
				eq_Prev[col,row] = limit(
					map(x->PrevN(x,H),sol.u)
				)
			end
		end
	end
end;

# ╔═╡ d1a88f76-1176-11ec-2a61-8d9ea1b23669
S = [[limit([sol.u[i][j] for i in 1:length(sol)]) for j in 1:2] for sol in Sol[1,:]]

# ╔═╡ b519decc-153a-11ec-1969-3548bc3f7c2e
begin
	scatter([S[i][1] for i ∈ 1:length(S)],[S[i][2] for i ∈ 1:length(S)],legend=false,xlabel="X",ylabel="Y")
end

# ╔═╡ cd36088a-ff54-11eb-286c-d580e2400678
begin
	if plot_eq
		p = plot(range(0.0,stop=2*N,length=100),[eq_Load[N,:] .- eq_Load_ill[N,:], eq_Load_ill[N,:]],label=["Eqilibrium Mutation Load (healthy)" "Equilibrium Mutation Load (ill)"], legend=:topleft, xlabel="De novo influx",framestyle = :zerolines,title="N=$N")
		
#		if N <= 2
#			plot!(p,
#				csvN[N].dni,
#				[csvN[N].mut,csvN[N].prev],
#				label=["Mutation Load" "Prevalence"]
#			)
#		end
		savefig(p,"PlotOutput/LoadIll_Eq_N=$N.pdf")
		p
	end
end

# ╔═╡ 7ebab978-0029-11ec-3df4-830eaeff284e
let
	if plot_eq
		p = plot(legend=:bottomright, 
			xlabel="De novo influx",
			framestyle = :zerolines,
			title="Eqilibrium Mutation Load"
		)

		for n in Ns

			plot!(p,
				range(0.0,stop=2*n,length=100),
				eq_Load[n,:],
				label="N=$n", 
			)
		end

		line = 0.0:0.1:2*Ns[end]
		plot!(p,line,line,linestyle=:dash,color=:black,label="")

		savefig("PlotOutput/Eq_Mut_Load.pdf")

		p
	end
end

# ╔═╡ df9ca870-0029-11ec-1dfb-4f48880a6369
let
	if plot_eq
		p = plot(legend=:bottomright, 
			xlabel="De novo influx",
			framestyle = :zerolines,
			title="Eqilibrium Prevalence"
		)

		for n in Ns

			plot!(p,
				range(0.0,stop=2*n,length=100),
				eq_Prev[n,:],
				label="N=$n", 
			)
		end

		savefig("PlotOutput/Eq_Prev.pdf")

		p
	end
end

# ╔═╡ 9d10a26e-0413-11ec-2c72-cd5cc8eb0a0f
let
	cmax = 2.1 .^(- Ns .+ 1)
	
	anim = @animate for i ∈ 1:length(dnis)
	
		l = @layout [
			grid(1,length(Ns))
   	 		b{0.2h}
			]
		
		Hs = []
		
		for N ∈ Ns
			h = heatmap(
					reshape(reverse(Sol[N,i].u[end]),3^N,1),
					clim=(0,cmax[N]),
					axis = false,
					ticks = false,
					cbar = false
				)
			push!(Hs,h)
		end
		push!(Hs, plot_dni_bar(i))
		plot(Hs...,
			layout = l,
			size = (600,600),
			title=reshape(push!(["N=$N" for N ∈ Ns],""),1,length(Ns)+1)
			)
	end
	#gif(anim, "PlotOutput/DetailEvol_N=1,_,$(Ns[end]).gif", fps = 15)	
end

# ╔═╡ 118c8ac8-fe6e-11eb-3ca8-8f13b3870c0d
let
	
	H = healthy_config_positions(N)
	H̄ = 1 .- H
	W = Load_weights(N)
	
	plot(
		sol.t,[
			map(x->LoadN(x,W), sol.u),
			map(x->LoadN(x .* H,W), sol.u),
			map(x->LoadN(x .* H̄,W), sol.u)
			#map(x->PrevN(x,healthy_config_positions(N)), sol.u)
			],
		label=["Total Mutation Load" "Load of Healthy Subpopulation" "Load of Ill Subpopulation"
		]
	)
	#hline!([limit(Load.(sol.u)),limit(Prev.(sol.u))])
end


# ╔═╡ 4176c240-03fe-11ec-3652-41d3308d1ebb
let
	Ns = 1:10
	ps=[]
	
	for N in Ns
		p = heatmap(
				reshape(reverse(healthy_config_positions(N)),3^N,1),
				c=cgrad([:black, :lightgray]),
				axis=false,
				ticks=false,
				cbar=false
				)
		push!(ps,p)
	end
	
	plot(ps..., 
		layout = (1,Ns[end]), 
		title=reshape(["$N" for N ∈ Ns],1,length(Ns)),
		size=(2400,1800)
		)
	
	savefig("PlotOutput/HealthyConfigPositionN=1,_,10.png")
end

# ╔═╡ 4b064d90-0650-11ec-1348-69242b2b875c
let
	Ns = 1:20
	q(n) = 1-(2^n)/(3^n)
	plot(Ns,q.(Ns))
end

# ╔═╡ Cell order:
# ╠═a9dcf66c-00e1-11ec-2489-85a8bf50529f
# ╠═def9fb92-00e6-11ec-0ce0-8bf38dfba80b
# ╠═b1cb8c58-00e1-11ec-060f-7bc5f20723c6
# ╠═6811d1c2-fe6a-11eb-05f8-151a60d2bd0a
# ╠═e2278904-37e3-11ec-07e6-578d2e0e155c
# ╠═91604328-ff38-11eb-0441-595f7d6f456f
# ╠═d1a88f76-1176-11ec-2a61-8d9ea1b23669
# ╠═b519decc-153a-11ec-1969-3548bc3f7c2e
# ╟─e002f598-fcd6-11eb-0693-cd37b6731ddb
# ╟─e2864ca2-00eb-11ec-15af-13e2c905c2d9
# ╠═118c8ac8-fe6e-11eb-3ca8-8f13b3870c0d
# ╠═19aeceb6-0bcd-11ec-2237-a19898d1df26
# ╟─2dd78c3a-00e6-11ec-33a2-07189bcd1803
# ╟─edf4f43a-01c7-11ec-3f89-87267d95ba52
# ╟─45ee92e0-00ee-11ec-2060-c791823a5884
# ╟─478f641c-00ec-11ec-1f4a-cd703ab756d6
# ╟─634cae54-01c7-11ec-34df-9df3b4d9558d
# ╠═cd36088a-ff54-11eb-286c-d580e2400678
# ╟─7ebab978-0029-11ec-3df4-830eaeff284e
# ╟─df9ca870-0029-11ec-1dfb-4f48880a6369
# ╟─9d10a26e-0413-11ec-2c72-cd5cc8eb0a0f
# ╟─ff76b4b2-04e4-11ec-1e43-59b314d25a53
# ╠═c6b86160-04e5-11ec-050c-711567fe7ca8
# ╟─5ca4273e-01c8-11ec-078a-5798323bc39d
# ╟─0eb45f48-ff58-11eb-2704-4f5dad08241b
# ╟─c4fdb038-00f7-11ec-3011-c9183ec5fe36
# ╠═b1a05bb4-00f0-11ec-1a54-13b82c896a6c
# ╟─41c242fe-ff58-11eb-01b5-81fac3c399cd
# ╠═18afc95a-fcca-11eb-2c88-792999cb5897
# ╟─abb776ea-00cf-11ec-3d9a-dbb357f9883e
# ╟─a445432e-fce7-11eb-33f3-b71379b6eddb
# ╠═685cf83e-fce7-11eb-38ee-cb456fdbc26d
# ╠═f3d2bb1c-00d0-11ec-3fc6-b7db47f9232e
# ╟─acdef9b4-fdf3-11eb-0aaa-6de346bde1cd
# ╟─8e4d6f14-fdef-11eb-1e04-278b64849bff
# ╠═670bcc9c-01d0-11ec-16a2-417c2d8d1924
# ╟─7d72ff6e-03f4-11ec-1b98-8f66c650f673
# ╟─71df4b56-fe73-11eb-3559-55a35d29f37c
# ╠═73f9e942-01cc-11ec-3482-d58071c89d6e
# ╟─5e6f439c-fe8e-11eb-1e33-2f7be0eaadff
# ╠═95801cc2-ff37-11eb-28c2-b16f7958283a
# ╠═c0d9a06c-fe78-11eb-3776-1723db5d3551
# ╟─6b6402f6-fcde-11eb-0a18-63cc8c83fae8
# ╟─27cbec3a-01d3-11ec-01c1-a7495311e4af
# ╠═f7c16cac-fcde-11eb-1a1e-37cbc0e31289
# ╟─bedb15aa-fced-11eb-0f8c-6b0f1063b26f
# ╟─df0f60ec-fdf1-11eb-23e0-9f2747d9d223
# ╟─a01907ac-03ee-11ec-3f15-81f0a341a1bb
# ╠═4176c240-03fe-11ec-3652-41d3308d1ebb
# ╠═4b064d90-0650-11ec-1348-69242b2b875c
# ╠═38ad5806-f5c2-11eb-1a4b-757cb20b199c

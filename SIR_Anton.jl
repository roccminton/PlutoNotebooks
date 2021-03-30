### A Pluto.jl notebook ###
# v0.12.10

using Markdown
using InteractiveUtils

# ╔═╡ 9f18ccd0-8e25-11eb-20e6-ab941f663550
using Plots

# ╔═╡ 10a9374c-8e1a-11eb-10b4-fb8272dd46d1
md"""
## SIR Model Simulation

#### Intro

The **SIR Model** has been developed to simulate an epidemic over time. The model consists of a system of 3 differential equations that express the rates of change of 3 variables over time. The 3 variables are:

1. **S** - the susceptibles of getting the infection
2. **I** - the infected
3. **R** - the recovered from the infection
3. D - the number of dead

#### The Model
Here follow the 3 equations that govern the model dynamics:

${ \begin{align*}
\frac{\mathrm{d}S}{\mathrm{d}t} & = -\beta \cdot S \cdot I + \delta \cdot R \\
\frac{\mathrm{d}I}{\mathrm{d}t} & = \beta \cdot S \cdot I + \gamma \cdot I \\
\frac{\mathrm{d}R}{\mathrm{d}t} & = \gamma \cdot I - \delta \cdot R \\
\end{align*} }$

#### Model simulation
The following is a simulation of the model described above
"""

# ╔═╡ 43970fa0-8e1b-11eb-3052-1701966a4478
begin
	#Model parameters
	beta1 = 0.14 	#rate of infection
	gamma1 = 0.07 	#rate of recovery (try also 0.07)
	delta1 = 0.00 	#rate of immunity loss
	rho1 = 0.001 	#rate of death of infected
	N1 = 6*10^7 	#Total population N = S + I + R + D
	I0 = 100.0 	#initial number of infected
	beta2 = 0.0 	#rate of infection
	gamma2 = 0.035 	#rate of recovery (try also o.07)
	delta2 = 0.0 	#rate of immunity loss
	rho2 = 0.2 	#rate of death of infected
	beta3 = 0.0 	#rate of infection from group 1 to group 2
	beta4 = 0.0 	#rate of infection from group 2 to group 1
	N2 = 1*10^6 #Total population N = S + I + R + D
	J0 = 0.0 	#initial number of infected
	gamma = 0.0 	
	T = 38*30 	#period of 300 daysplot
	dt = 1/400 	#time interval of 6 hours (1/4 of a day)
	
	md"""
	Value of parameter R10 is $(beta1/gamma1) | Value of parameter R20 is $(beta2/gamma2)
	"""
end

# ╔═╡ 883a4ad8-8e26-11eb-06c4-e30f3c0f70c1
# Plots that display the epidemic outbreak
tt = 0:dt:T-dt;

# ╔═╡ 5f948160-8e27-11eb-1003-fbcf284800e6
#Map
#plot(tt,[I1,I2,D1,D2],color=[:blue :red :green :yellow], lw=2, grid=true, xlabel="Days", ylabel="Number of individuals", label=["I1" "I2" "D1" "D2"])

# ╔═╡ c30b3fb4-8e1c-11eb-05c7-d7f95aa97fcf
md"Function that calculates the model evolution at time period T"

# ╔═╡ 01a236b0-8e1d-11eb-2630-ebff3f351bb3
function sir_model(beta1,gamma1,delta1,rho1,N1,beta2,gamma2,delta2,rho2,N2,beta3,beta4,I0,J0,T,dt)
	S1 = zeros(Int(T/dt))
	S1[1] = N1
	I1 = zeros(Int(T/dt))
	I1[1] = I0
	R1 = zeros(Int(T/dt))
	D1 = zeros(Int(T/dt))
	S2 = zeros(Int(T/dt))
	S2[1] = N2
	I2 = zeros(Int(T/dt))
	I2[1] = J0
	R2 = zeros(Int(T/dt))
	D2 = zeros(Int(T/dt))
	I = zeros(Int(T/dt))
	J1 = zeros(Int(T/dt))
	J2 = zeros(Int(T/dt))
	gamma = gamma1
	
	for tt ∈ 1:Int(T/dt)-1
		if I2[tt] + I1[tt] > 10
			if tt > 10/dt
				if J2[tt-Int(5/dt)] > 200
					gamma = 2.8*gamma1
				end
				if J2[tt-Int(5/dt)] < 17
					gamma = gamma1
				end
			end

			dS1 = (-beta1*I1[tt]*S1[tt]/N1+delta1*R1[tt]-beta4*S1[tt]*I2[tt]/N1) * dt
			dI1 = (beta1*I1[tt]*S1[tt]/N1-gamma*I1[tt]+beta4*S1[tt]*I2[tt]/N1) * dt
			dR1 = (gamma*(1-rho1)*I1[tt]-delta1*R1[tt]) * dt
			dD1 = (gamma*rho1*I1[tt]) * dt
			dJ1 = beta1*I1[tt]*S1[tt]/N1

			S1[tt+1] = S1[tt] + dS1
			I1[tt+1] = I1[tt] + dI1
			R1[tt+1] = R1[tt] + dR1
			D1[tt+1] = D1[tt] + dD1
			J1[tt+1] = J1[tt] + dJ1

			dS2 = (-beta2*I2[tt]*S2[tt]/N2+delta2*R2[tt]-beta3*I1[tt]*S2[tt]/N1) * dt
			dI2 = (beta2*I2[tt]*S2[tt]/N2-gamma2*I2[tt]+beta3*I1[tt]*S2[tt]/N1) * dt
			dR2 = (gamma2*(1-rho2)*I2[tt] - delta2*R2[tt]) * dt
			dD2 = (gamma2*rho2*I2[tt]) * dt

			S2[tt+1] = S2[tt] + dS2
			I2[tt+1] = I2[tt] + dI2
			R2[tt+1] = R2[tt] + dR2
			D2[tt+1] = D2[tt] + dD2
			if tt > 7/dt +1
				#I[tt+1] = I1[tt+1]+I2[tt+1]+R1[tt+1]+R2[tt+1]-(I1[tt-7/dt]+I2[tt-7/dt]+R1[tt-7/dt]+R2[tt-7/dt])
				#J1[tt+1] = (I1[tt+1]+R1[tt+1]+D1[tt+1]-D1[tt-7/dt]-(I1[tt-7/dt]+R1[tt-7/dt]))*100000/N1
				J2[tt+1] = (J1[tt+1] - J1[tt-Int(7/dt)])*100000/N1
				#J2[tt+1] = I2[tt+1]+R2[tt+1]+D2[tt+1]-D2[tt-7/dt]-(I2[tt-7/dt]+R2[tt-7/dt])
			end
		else
			print("boom")
			break
		end
	end
	return S1, I1, R1, D1, S2, I2, R2, D2, I, J1, J2
end;

# ╔═╡ 81f4c3f4-8e23-11eb-2f4e-f1b0173d3ef6
#Calculate the model
S1, I1, R1, D1, S2, I2, R2, D2, I, J1, J2 = sir_model(beta1,gamma1,delta1,rho1,N1,beta2,gamma2,delta2,rho2,N2,beta3,beta4,I0,J0,T,dt);

# ╔═╡ 116cc80c-8e25-11eb-1120-911f0eae6a3c
# Curves
plot(tt,[S1,I1,R1,D1],color=[:blue :red :green :yellow], label=["S" "I" "R" "D"], lw=5, grid=true,xlabel="Days",ylabel="Number of individuals")

# ╔═╡ c201ce30-8e26-11eb-15a1-99f683ff6d20
plot(tt,[S2,I2,R2,D2],color=[:blue :red :green :yellow], label=["S2" "I2" "R2" "D2"], lw=2, grid=true,xlabel="Days",ylabel="Number of individuals")

# ╔═╡ f7b43f62-8e27-11eb-1675-55e0b4f3462e
#plot incidence
plot(tt,J1,color=:blue,lw=2,xlabel="Days",ylabel="7 day incidence",label="I")

# ╔═╡ 5f5686ac-8e28-11eb-24a1-299e86f64808
#plot incidence
plot(tt,J2,color=:blue,lw=2,xlabel="Days",ylabel="7 day incidence",label="I")

# ╔═╡ Cell order:
# ╟─10a9374c-8e1a-11eb-10b4-fb8272dd46d1
# ╠═43970fa0-8e1b-11eb-3052-1701966a4478
# ╠═81f4c3f4-8e23-11eb-2f4e-f1b0173d3ef6
# ╠═883a4ad8-8e26-11eb-06c4-e30f3c0f70c1
# ╠═116cc80c-8e25-11eb-1120-911f0eae6a3c
# ╠═c201ce30-8e26-11eb-15a1-99f683ff6d20
# ╠═5f948160-8e27-11eb-1003-fbcf284800e6
# ╠═f7b43f62-8e27-11eb-1675-55e0b4f3462e
# ╠═5f5686ac-8e28-11eb-24a1-299e86f64808
# ╟─c30b3fb4-8e1c-11eb-05c7-d7f95aa97fcf
# ╠═01a236b0-8e1d-11eb-2630-ebff3f351bb3
# ╠═9f18ccd0-8e25-11eb-20e6-ab941f663550

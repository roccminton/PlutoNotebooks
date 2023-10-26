### A Pluto.jl notebook ###
# v0.19.26

using Markdown
using InteractiveUtils

# ╔═╡ eeab1982-e223-40bd-a95b-ac1178fdaea7
md"## Blatt 3"

# ╔═╡ 5485d046-9521-4f4b-8a16-1294a219433c
md"#### Aufgabe 1"

# ╔═╡ d209d312-7423-11ee-26e2-833b87105abb
begin
	#Helping Function returning the sum of ones for a given type
	function two(::Type{T}) where T<:AbstractFloat
		return one(T) + one(T)
	end
	#General Function for calculating the machine epsilon for a given input type
	function my_eps(::Type{T}) where T<:AbstractFloat
		ε = one(T)
		ε_old = one(T)
		while one(T) + ε ≠ one(T)
			ε_old = ε
			ε = T(ε_old) / two(T)
		end
		return ε_old
	end
	#Default version calculating the machine epsilon of a Float64
	my_eps() = my_eps(Float64)
end;

# ╔═╡ d0fadf2e-398e-4d2e-848a-8fa1a5b31520
md""" My machine epsilon calculatetd in the above way with the `my_eps()` function results in 
"""

# ╔═╡ 46ed26c3-4db7-4d21-bf66-cf17800b731b
 my_eps()

# ╔═╡ be087d6f-6d65-4893-8bff-1ec4294ca6f3
md"whereas the machine epsilon calculated by Julia with the function `eps()` results in" 

# ╔═╡ 633fa02c-8d72-4e61-bc25-edcc298983f7
eps()

# ╔═╡ 2f1cff9c-aa67-4a6d-b6ad-adf46cdaec8d
md"Pretty cool! Seems that Julia uses the same method to calculate the machine epsilon 👍"

# ╔═╡ 447a9937-9459-4c36-b153-36bb54c93029
md"---"

# ╔═╡ 289b8f33-93a4-4bd7-ba27-78209de35a74
md"#### Aufgabe 2"

# ╔═╡ 688a5a1f-852d-4627-99b8-e031eed24f8f
md"First define some helper functions for a better looking notebook. These are not necessary but pretty"

# ╔═╡ 59a16ca9-5959-4011-b309-9271d8c1f322
begin
	#Function that converts a 8-bit Array into a decimal
	eightbit_to_dec(d) = (-1)^d[1] * sum(d[8-i]*2.0^(i-3) for i in 0:6)	
	#Function that converts a 8-bit Array to a string
	eightbit_to_string(d) = string(round.(Integer,d[1:5])...,".",round.(Integer,d[6:8])...)
	#Function that converts a string to 8-bit Array
	string_to_eightbit(s) = [parse(Int,c,base=10) for c in filter(x->x!='.',s)]
end;

# ╔═╡ eb858b52-36d2-4f65-8f88-7608991c1de1
md"From these we immediately get the $x_{min} > 0$ and $x_{max} > 0$ since the biggest positive number that can be represented within this system is $01111.111$ and the smallest positive number is $00000.001$"

# ╔═╡ 30fbcf73-6288-4fb0-b8e8-1c20301b0687
x_min = eightbit_to_dec([0,0,0,0,0,0,0,1])

# ╔═╡ eb87aed6-1b26-40d9-8fa6-724fb6cba3c2
x_max = eightbit_to_dec([0,1,1,1,1,1,1,1])

# ╔═╡ 76a464d1-86f2-4418-8bcc-9753b36cd957
md"Then define a function that can convert your numbers into arrays of size 8"

# ╔═╡ db24cdc7-5134-4201-833f-aaafcad8e792
#Function that converts a number to a eightbith vector and returns both the vector and the remainder if the number cannot be split evently
function number_to_eightbit(a)
	d = zeros(Integer,8)
	a < 0 && (d[1] = 1 ; a *= (-1))
	a > x_max && error("The input $a is too big")
	for i in 6:-1:0
		if a-2.0^(i-3) < 0 
			d[8-i] = 0
		else
			d[8-i] = 1
			a -= 2.0^(i-3)
		end
	end

	return d,a
end

# ╔═╡ 9e9b9975-338c-4110-8be5-7755908c6a5d
md"And that all we need to convert the numbers from part a)"

# ╔═╡ dc48e2b6-c353-4d2c-98d5-93819b1f6a0e
Markdown.parse("7.25 = " * eightbit_to_string(number_to_eightbit(7.25)[1]))

# ╔═╡ c43620b7-43fb-4844-ae29-6f5a65b6935b
Markdown.parse("-5.625 = " * eightbit_to_string(number_to_eightbit(-5.625)[1]))

# ╔═╡ c4e68ff8-0956-47e5-98cb-2b469f61884f
md"Obiously one coud have done this by hand, which would have been much faster... **but** this way it was much more fun, you don't have to compute nothing in your head and you can easily convert any number (given its in the range, which brings us to part c)"

# ╔═╡ 577bdf2d-f061-4901-88a8-6ae50e779502
x = 1/3

# ╔═╡ 4b082573-b0c4-4cd0-a22e-04f2cd2e9208
md"Fortunately we already hand out the absolute errow as the remainder of the calculation. We can see that they are equal if we convert the number once in the dual and then back. Let call this number y."

# ╔═╡ a50fa85f-e38e-42f5-92a2-f9413a1447e7
y = eightbit_to_dec(number_to_eightbit(x)[1])

# ╔═╡ 6c13fe6b-c26c-4911-b903-ec0fb3e6ca7c
md"Then we get the absolute error as"

# ╔═╡ f6cadca2-e98d-49e7-98b7-d70a41409a9a
error_abs = abs(x-y)

# ╔═╡ 37ad8f61-9b62-4bb5-881e-f4c617d836ea
md"which is the same as the remainder of the calculation above"

# ╔═╡ 85dc7e46-aa03-4b68-8081-102d6f03ad82
number_to_eightbit(x)[2]

# ╔═╡ 08d66d48-b98f-402b-b960-e0ed194f2c8d
md"Let's put this in a function as well, maybe we can use it later on"

# ╔═╡ 79ec646e-baf4-4495-afb0-dc9d256c543f
ε_abs(x) = number_to_eightbit(x)[2]

# ╔═╡ e1c0cad9-c54b-48fe-be7c-f8ec44a8cc6b
md"Moreover the relative error is given by"

# ╔═╡ 0e897861-7d5b-43f8-ab1f-aa65e97e8bcb
error_rel = error_abs / abs(y)

# ╔═╡ 7c02051d-07f3-4d61-941c-704a168d6f57
md"Let's put also this calculation in a function... you'll never know."

# ╔═╡ e95b0b12-d17a-474d-bd5c-6f48bc9d6baf
ε_rel(x) = ε_abs(x) / abs(eightbit_to_dec(number_to_eightbit(x)[1]))

# ╔═╡ 307aac5f-3db4-4b2d-8cf4-efa5ce07c642
md"The maximal absolute error occurs when we are juuuuust above the minimal step size of $2^{-3}$ that we can within this system. Luckily we just calculated the most small number we can represent with this machine in exercise 1. So lets see what the maximal absolute error is"

# ╔═╡ 6941f869-da48-4be9-a017-db76dd09c0e3
ε_abs_max = ε_abs(2^(-3)-my_eps())

# ╔═╡ 00bed98f-d246-4789-a74a-c416426489ba
md"The relativ error on the other hand converges to one as we get closer and closer to $2^{-2}$ from below as one can see here"

# ╔═╡ 261619a7-0361-4b3e-b656-01386678e927
ε_rel_max = ε_rel(2^(-2)-my_eps())

# ╔═╡ 79d00237-f51a-41fe-b3e2-02e08488f5e4
md"---"

# ╔═╡ Cell order:
# ╟─eeab1982-e223-40bd-a95b-ac1178fdaea7
# ╟─5485d046-9521-4f4b-8a16-1294a219433c
# ╠═d209d312-7423-11ee-26e2-833b87105abb
# ╟─d0fadf2e-398e-4d2e-848a-8fa1a5b31520
# ╠═46ed26c3-4db7-4d21-bf66-cf17800b731b
# ╟─be087d6f-6d65-4893-8bff-1ec4294ca6f3
# ╠═633fa02c-8d72-4e61-bc25-edcc298983f7
# ╟─2f1cff9c-aa67-4a6d-b6ad-adf46cdaec8d
# ╟─447a9937-9459-4c36-b153-36bb54c93029
# ╟─289b8f33-93a4-4bd7-ba27-78209de35a74
# ╟─688a5a1f-852d-4627-99b8-e031eed24f8f
# ╠═59a16ca9-5959-4011-b309-9271d8c1f322
# ╟─eb858b52-36d2-4f65-8f88-7608991c1de1
# ╠═30fbcf73-6288-4fb0-b8e8-1c20301b0687
# ╠═eb87aed6-1b26-40d9-8fa6-724fb6cba3c2
# ╟─76a464d1-86f2-4418-8bcc-9753b36cd957
# ╠═db24cdc7-5134-4201-833f-aaafcad8e792
# ╟─9e9b9975-338c-4110-8be5-7755908c6a5d
# ╟─dc48e2b6-c353-4d2c-98d5-93819b1f6a0e
# ╟─c43620b7-43fb-4844-ae29-6f5a65b6935b
# ╟─c4e68ff8-0956-47e5-98cb-2b469f61884f
# ╠═577bdf2d-f061-4901-88a8-6ae50e779502
# ╟─4b082573-b0c4-4cd0-a22e-04f2cd2e9208
# ╠═a50fa85f-e38e-42f5-92a2-f9413a1447e7
# ╟─6c13fe6b-c26c-4911-b903-ec0fb3e6ca7c
# ╠═f6cadca2-e98d-49e7-98b7-d70a41409a9a
# ╟─37ad8f61-9b62-4bb5-881e-f4c617d836ea
# ╠═85dc7e46-aa03-4b68-8081-102d6f03ad82
# ╟─08d66d48-b98f-402b-b960-e0ed194f2c8d
# ╠═79ec646e-baf4-4495-afb0-dc9d256c543f
# ╟─e1c0cad9-c54b-48fe-be7c-f8ec44a8cc6b
# ╠═0e897861-7d5b-43f8-ab1f-aa65e97e8bcb
# ╟─7c02051d-07f3-4d61-941c-704a168d6f57
# ╠═e95b0b12-d17a-474d-bd5c-6f48bc9d6baf
# ╟─307aac5f-3db4-4b2d-8cf4-efa5ce07c642
# ╠═6941f869-da48-4be9-a017-db76dd09c0e3
# ╟─00bed98f-d246-4789-a74a-c416426489ba
# ╠═261619a7-0361-4b3e-b656-01386678e927
# ╟─79d00237-f51a-41fe-b3e2-02e08488f5e4

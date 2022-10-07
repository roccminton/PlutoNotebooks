### A Pluto.jl notebook ###
# v0.17.1

using Markdown
using InteractiveUtils

# ╔═╡ 463cb478-327b-431f-a2f9-af2da22d1a24
md"#### Sudoku-Sums"

# ╔═╡ f1678f15-2c97-4048-b544-0236622f2388
md"""

The function $\text{sums(res,digits)}$ plots all the possible ways to produce the result $\text{res}=s_1 + \cdots +s_n$ as a sum of $n=\text{digits}$ where the summands are digits from one to nine $s_1,\cdots,s_n \in \{1,\cdots,9\}$ and cannot appear more than once $s_i \neq s_j$ for $i \neq j$.

One can pass additional arguments to specify which digits have to be contained as summand or alternatively which cannot be contained as summand. In that case one has to add another additional argument $\text{not=true}$.

"""

# ╔═╡ b0d9e534-fda6-425d-8e37-74ba246c0ee5
function emptyintersect(S,T)
	P = Set()

	for s in S
		for t in T
			isempty(s∩t) && push!(P,(s,t))
		end
	end
	return P
end

# ╔═╡ f9b881bb-ea63-4424-8872-8383f6412f29
begin
	function maxsum(ndigits)
		if 1 <= ndigits <= 9
			return sum(i for i in 10-ndigits:9)
		else
			return 0
		end
	end

	function minsum(ndigits)
		if 1 <= ndigits <= 9
			return sum(i for i in 1:ndigits)
		else
			return 0
		end
	end
end

# ╔═╡ 9fd7d9b9-f1c9-4a97-8241-7835fdb7f5aa
begin
	function sums(res,digits)
		digits==2 && return sums(res)
		ways = Set()
		(res < minsum(digits) || res > maxsum(digits)) && return ways
		(digits < 1 || digits > 9) && return ways

		min, max = minsum(digits-1),maxsum(digits-1)
		for i in 1:9
			if min ≤ res-i ≤ max
				way = filter(S->length(S)==digits,push!.(sums(res-i,digits-1),i))
				!isempty(way) && push!(ways,way...)
			end
		end
		return ways
	end

	function sums(res)
		ways = Set()
		(res < 3 || res > 17) && return ways
		for i in 1:9, j in 1:9
			if i < j && i+j == res
				push!(ways,Set([i,j]))
			end
		end
	return ways
	end

	sums(res,digits,c::Real,not=false) = filter!(
		W->(not ? (c∉W) : (c∈W)),
		sums(res,digits)
	)
	sums(res,digits,c::Set,not=false) = filter!(
		not ? W->all(cs∉W for cs∈c) : (W->c⊆W),
		sums(res,digits)
	)	
	sums(res,digits,c::Vector,not=false) = sums(res,digits,Set(c),not)
end

# ╔═╡ Cell order:
# ╟─463cb478-327b-431f-a2f9-af2da22d1a24
# ╟─f1678f15-2c97-4048-b544-0236622f2388
# ╠═b0d9e534-fda6-425d-8e37-74ba246c0ee5
# ╠═9fd7d9b9-f1c9-4a97-8241-7835fdb7f5aa
# ╠═f9b881bb-ea63-4424-8872-8383f6412f29

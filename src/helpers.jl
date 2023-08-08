import Base.zeros, Base.println
println(cond::Bool, io::IO, xs...) = if cond; println(io, xs...) end
println(cond::Bool, xs...)         = if cond; println(xs...) end

zeros(x::Float64)            = zeros(Int(round(x)))
zeros(x::Float64,y::Int64)   = zeros(Int(round(x)), y)
zeros(x::Float64,y::Float64) = zeros(Int(round(x)), Int(round(y)))

vecF64(x::Any) = Vector{Float64}(vec(x))
vecI64(x::Any) = Vector{Int64}(vec(x))
vecC64(x::Any) = Vector{ComplexF64}(vec(x))

nan_to_zero(v) = map(x -> isnan(x) ? zero(x) : x, v)
nan_to_inf(v)  = map(x -> isnan(x) ? -Inf : x, v)

function smooth!(v::Matrix)
	for i=1:size(v,1), j=1:size(v,2)
		if isnan(v[i,j])
			v[i,j] = if j == 1
				v[i,j+1]
			elseif j == size(v,2)
				v[i,j-1]
			else
				(v[i,j-1] + v[i,j+1]) / 2.0
			end
		end
	end
	return v
end

function vals(d::Dict{Symbol,Float64})
	return [d[x] for x in keys(d)]::Vector{Float64}
end
vals(d::OrderedDict{Symbol,T}) where T<:Number = [d[x] for x in keys(d)]::Vector{T}

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T<:Float64}
  m, n = length(v1), length(v2)
  v1   = reshape(v1, m, 1)
  v2   = reshape(v2, 1, n)
  return repeat(v1, 1, n), repeat(v2, m, 1)
end

"""
```
build_θ(θ_free_val::Vector{S}, θ_fix_val::Vector{S}, free_ind::Vector{T},
		fix_ind::Vector{T}) where {S<:Float64, T<:Int64}
```
This reconstitutes the full parameter vector by slotting the free and fixed
parameters back into their original indices.

# Motivation
The reason we extract the fixed parameters from the vector pre-optimization is
for performance: depending on the optimization algorithm, if the optimizer
erroneously believes the space to optimize over is +`N_fix` dimensions larger,
efficiency of exploration may sharply decline.
"""
function build_θ(θ_free_v::Vector{S}, θ_fix_v::Vector{S}, free_i::Vector{T},
				 fix_i::Vector{T}) where {S<:Float64, T<:Int64}
	N_θ = length(fix_i) + length(free_i)
	θ   = zeros(N_θ)
	θ[free_i] .= θ_free_v
	θ[fix_i]  .= θ_fix_v
	return θ
end

"""
```
demand_shopper(α::T, β::T, p::V, cdid::V, obs_w::V;
               testing::Bool=false) where {T<:Float64, V<:Vector{Float64}}
```
Solve for the demand of shoppers, along with 1st & 2nd order derivatives.

Note utility for outside good -> price that is the lowest (leave condition out).
"""
function demand_shopper(α::T, β::T, p::V, cdid::Vector{Int64}, obs_w::V;
                        α_c::T = 0.0, cond::V = V(undef, length(cdid)),
						min_p::T = 0.0, has_min_p::V = V(undef, length(cdid)),
						numlist::V = V(undef, length(cdid)),
						allout::Bool = false) where {T<:Float64, V<:Vector{T}}
	N    = length(cdid)
    expU = @. exp((p * -α) + (cond * -α_c) + (min_p * has_min_p))
    sum1 = sum(sparse(collect(1:N), cdid, obs_w .* expU); dims=1)' .+ exp(β * α)
    sum2 = sum1[cdid]

    D0   = @.       expU                                     /  sum2
    dD0  = @. -α  * expU * (sum2 - expU)                     / (sum2 ^ 2)
    d2D0 = @. α^2 * expU * (sum2 - expU) * (sum2 - (expU*2)) / (sum2 ^ 3)

    replace!.([D0, dD0, d2D0], NaN => 0)

    if allout
        return D0, dD0, d2D0, sum2, expU
    else
        return D0, dD0, d2D0
    end
end

"""
```
solve_γ(p_in::V, D0::V, dD0::V, d2D0::V, δ::T, η::T, γ0::V, r::T, ϵ::T;
        allout = false) where {T<:Float64, V<:Vector{Float64}}
```
Wrapper for solving for the γ rationalizing price choice.
"""
function solve_γ(p::Union{U,V}, D0::Union{U,V}, dD0::Union{U,V}, d2D0::Union{U,V},
                 δ::T, η::T, γ0::Union{U,V}, r::T, ϵ::T;
				 η_c::T = 0.0, cond::V = V(undef, length(cdid)),
                 allout=false) where {T<:Float64, V<:Vector{T}, U<:Vector{ComplexF64}}

    Dm   = @. δ * (cond ^ (-η_c)) * (p ^ (-η))
    dDm  = @. δ * (cond ^ (-η_c)) * (p ^ (-η-1)) * (-η)
    d2Dm = @. δ * (cond ^ (-η_c)) * (p ^ (-η-2)) * (η) * (η+1)

    # Adjust for rounding error
    D0_ϵ  = @. D0  + ϵ * dD0 + 0.5 * d2D0 * ϵ^2
    dD0_ϵ = @. dD0 + ϵ * d2D0

    if allout
        return solve_γ(Dm, D0_ϵ, dDm, dD0_ϵ, γ0, p, r), Dm, dDm, d2Dm
    else
        return solve_γ(Dm, D0_ϵ, dDm, dD0_ϵ, γ0, p, r)
    end
end

"""
```
solve_γ(Dm::V, D0::Union{V,T}, dDm::V, dD0::Union{V,T}, γ0::Union{V,T}, p::Vector{T},
        r::T; allout::Bool = false) where {T<:Float64, V<:Vector{Float64}}
```
Solve for Poisson rate γ rationalizing a store i's title k price choice. [p 37]
"""
function solve_γ(Dm::Union{T,U,V}, D0::Union{T,U,V}, dDm::Union{T,U,V},
                 dD0::Union{T,U,V}, γ0::Union{T,U,V}, p::Union{T,U,V}, r::T;
                 allout = false) where {T<:Float64, V<:Vector{T}, U<:Vector{ComplexF64}}

    A  = @. Dm ^ 2.0
    B  = @. r * ((p * dDm) + Dm) + (2.0 * Dm * γ0 * D0)
    C  = @. (r * γ0) * ((p * dD0) +  D0) + (γ0 ^ 2) * (D0 ^ 2)

    γ1 = @. (-B + ((B ^ 2 - 4 * A * C) + 0im) ^ (0.5)) / (2 * A)

    if !allout; return γ1 end

    return γ1, (real.(γ1) .> 0) .& (imag.(γ1) .≈ 0)
end

"""
```
obs_cal(βσ3, data, basellh, demandcal, p0, ϵ, WFcal)
```
Corresponds to Masao/obscalnewtest2015.m. Works!

Start: (WFcal = false)
0.005568 seconds (37.75 k allocations: 3.399 MiB)
0.005290 seconds (33.48 k allocations: 2.851 MiB)
0.005211 seconds (33.48 k allocations: 2.670 MiB) -> Cut out extra returned

hcat:
56.179319 seconds (216.44 M allocations: 20.706 GiB, 30.97% gc time, 1.14% compilation time)
48.379337 seconds (216.41 M allocations: 19.084 GiB, 26.38% gc time, 1.38% compilation time)
39.416123 seconds (216.40 M allocations: 17.101 GiB, 14.07% gc time, 1.74% compilation time)
"""
function obs_cal(βσ3::V, #d::Dict{Symbol,Vector{<:Number}},
				 d_sym::Symbol, N::S, M::S, basellh::V, ϵ::T;
				 α_c::T = 0.0, η_c::T = 0.0, min_p::T = 0.0,
  				 demandcal::Bool = false, disap::V = Vector{Float64}(),
                 WFcal::Bool = false)::Vector{T} where {S<:Int64, T<:Float64,
                                                        V<:Vector{T}}
    # Conjure dictionary
	d::Dict{Symbol,Vector{<:Number}} = (d_sym == :d_on_09) ? d_on_09 :
	                                   (d_sym == :d_on_12) ? d_on_12 : bp
	@unpack numlist, d_first, p = d
	pdif = haskey(d, :pdif) ? d[:pdif] : d[:p]

	# Unpack and transform parameters
    γ0     = @. βσ3[1] * (numlist ^ βσ3[8]) / mean(numlist .^ βσ3[8])
    α      = (βσ3[2] + 1) * βσ3[14] ^ βσ3[15]
    β      = βσ3[3] / βσ3[14] ^ βσ3[15]
    m      = βσ3[4]
    η      = βσ3[6] + 1
    r      = βσ3[7]
    βcond  = βσ3[10]
    olp    = βσ3[13]
    δ      = βσ3[14]
    γscale = @. βσ3[5] / m * (numlist ^ βσ3[9] / mean(numlist .^ βσ3[9])) *
                exp(βσ3[12] * d[:localint])
    nat_disap = βσ3[16]

    # Solve for demand + its 1st & 2nd order derivatives
    D0, dD0, d2D0 = demand_shopper(α, β, pdif .- βcond .* d[:cond_dif] ./ α,
								   d[:cdid], d[:obs_w];
								   α_c = α_c,     cond = d[:condition],
								   min_p = min_p, has_min_p = d[:has_min_p])

    # Solve for (lower) γ rationalizing price choice
    p1 = p .- ϵ
    γ1 = solve_γ(p1, D0, dD0, d2D0, δ, η, γ0, r, -ϵ;
				 η_c = η_c, cond = d[:condition])


    # Solve for (upper) γ rationalizing price choice
    p2 = p .+ ϵ
    γ2, Dm, dDm, d2Dm = solve_γ(p2, D0, dD0, d2D0, δ, η, γ0, r, ϵ; allout=true,
								η_c = η_c, cond = d[:condition])

    SOC = @. r * p2 * (γ2 * d2Dm + γ0 * d2D0) + 2 * (r + (γ2 * Dm) +
             γ0 * (D0 + (ϵ * dD0) + 0.5 * (d2D0 * ϵ ^ 2))) *
            (γ2 * dDm + γ0 * (dD0 + ϵ * d2D0))

    γ1 = real.(γ1)
    γ2[imag.(γ2) .!= 0] .= 0
    γ2[real.(SOC) .> 0] .= 0
    γ2[real.(γ2)  .< 0] .= 0
    γ2 = real.(γ2)

    Π_2 = @. p2 / (r / (γ2 * Dm + γ0 * D0) + 1)
    Π_H = @. ((r * (η-1) / δ) ^ (-1/η) / (1/(η-1)+1)) * γ2 ^ (1/η)

    γ3 = solve_γ(Dm, 0.0, dDm, 0.0, 0.0, p2, r)

    γ2[real.(Π_H) .> real.(Π_2)] .= γ3[real.(Π_H) .> real.(Π_2)]
    γ1[real.(γ1) .< 0] .= .0
    γ1[real.(γ1) .> real.(γ2)] .= γ2[real.(γ1) .> real.(γ2)]

    γ1 = real.(γ1)
    γ2 = real.(γ2)

    #Dm = @. δ * (p ^ (-η))

    γ_dist = Gamma.(m, γscale)

    if demandcal == 1
        demandlh   = @. (disap > 0) - (2 * disap - 1) *
                        exp(-0.166666667 * (γ0 * D0 +
						(γ2 + γ1) / 2 * Dm)) * nat_disap
        demandlhol = @. (disap > 0) - (2 * disap - 1) *
                        # Probability of nondisappear due to shopper demand
                        exp(-0.166666667 * (γ0 * D0)) *
                        # Nonshopper demand (expectation taken wrt to γi)
                        (1 + γscale * 0.166666667 * Dm) ^ -m * nat_disap

        lip_o  = @. min.([cdf(γ_dist[i], γ2[i]) -
					      cdf(γ_dist[i], γ1[i]) for i=1:N], 1) * demandlh ^ 3
        lip_ol = @. basellh * demandlhol ^ 3 # Price likelihood
    else
        lip_o  = min.([cdf(γ_dist[i], γ2[i]) -
					   cdf(γ_dist[i], γ1[i]) for i=1:N], 1)
        lip_ol = basellh # Price likelihood
    end

    liptemp = @. (1 - olp) * lip_o + olp * lip_ol + 1e-10
    olppost = vec(olp .* lip_ol ./ liptemp)

	liptemp[liptemp .< 0.0] .= 0.0
    lip = log.(liptemp)

	if WFcal
		pi_v, CSns, CSs = welfaresimple(γ1, γ2, γscale .* m, γ0, olppost, Dm, D0,
                                 		pdif, p, N, M, d[:cdindex], d_first,
								 		[α; β; η; r])
    	#return [lip; γ2; γ1; γ0; D0; Dm; pi_v; CSns; CSs; olppost]
		return [lip; CSns; CSs; olppost; pi_v]
	else
		return lip
	end
end

"""
```
welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, pdif, data, scalarparas)
```
Corresponds to /Masao/welfaresimple.m. Works!
"""
function welfaresimple(γ1::V, γ2::V, γscale::V, γ0::V, olppost::V, Dm::V, D0::V,
                       pdif::V, p::V, N::S, M::S, cdindex::Vector{S},
                       d_first::Vector{S}, scalarparas::V;
					   N_draw = 10000) where {S<:Int64, V<:Vector{Float64}}
    α = scalarparas[1]
    β = scalarparas[2]
    η = scalarparas[3]
    r = scalarparas[4]

    γ1ave = 0.5 * (γ1 + γ2)

    pi_o  = @. p * (γ1ave  * Dm + γ0 * D0) / (r + γ1ave  * Dm + γ0 * D0)
    pi_ol = @. p * (γscale * Dm + γ0 * D0) / (r + γscale * Dm + γ0 * D0)
    pi_v  = @. olppost * pi_ol + (1 - olppost) * pi_o

    CSns_o  = @. (1/(η-1)) * γ1ave  * Dm * p / (r + γ1ave  * Dm + γ0 * D0)
    CSns_ol = @. (1/(η-1)) * γscale * Dm * p / (r + γscale * Dm + γ0 * D0)
    CSns    = @. olppost * CSns_ol + (1 - olppost) * CSns_o
    CSgain  = zeros(N)

    mktsize = @. cdindex - d_first + 1

    # for k = 1:M
    #     ind_k    = d_first[k]:cdindex[k]
	#
    #     best_p   = Vector{Float64}(undef, N_draw)
    #     best_ind = Vector{Float64}(undef, N_draw)
    #     β_p      = Vector{Float64}(undef, N_draw)
	#
    #     gumb_draw  = rand(Gumbel(0,1), mktsize[k] + 1, N_draw)
	#
    #     rand_price = repeat(-[pdif[ind_k]; -β], 1, N_draw) .- gumb_draw ./ α
	#
    #     best, bestindex = vec.(findmax(rand_price, dims = 1))
	#
    #     temp = sparse([(x->x[1]).(bestindex); mktsize[k]+1], 1:N_draw+1,
    #                   [best .- rand_price[end,:]; 1])
	#
    #     CSgain[ind_k] = sum(temp[1:mktsize[k], 1:N_draw],      dims=2) ./
    #                    (sum(temp[1:mktsize[k], 1:N_draw] .> 0, dims=2) .+ 1e-5)
    # end

    CSs_o  = @. (p/(η-1)) * γ0 * D0 / (r + γ1ave  * Dm + γ0 * D0)
    CSs_ol = @. (p/(η-1)) * γ0 * D0 / (r + γscale * Dm + γ0 * D0)
    CSs    = @. olppost * CSs_ol + (1 - olppost) * CSs_o

    return pi_v, CSns, CSs
end

"""
```
output_statistics(boot_out = "$OUTPUT/bootstrap_welfare.csv")
```
Outputs variables in same order of paper table.
"""
function output_statistics(; boot = Matrix{F64}(), distpara = Matrix{F64}(),
							 boot_out = "$OUTPUT/bootstrap_welfare.csv",
							 vint = "", write_out = false)

    if isempty(boot)
		boot = CSV.read(boot_out, DataFrame, header = false)
    	boot = unique(boot)
    	boot = boot[:,2:15]

		distpara = CSV.read(boot_out, DataFrame, header = false)
	    distpara = unique(distpara)
	    distpara = distpara[:,16:21]
	end

    N_bs   = size(boot, 1)
    b_boot = zeros(N_bs, 25)

    for i = 1:N_bs
        xx  =  Vector(boot[i, :])		# 1:14 -> boot_out[:, 2:15]
        est = [Vector(distpara[i, :]);	# 1:6  -> boot_out[:, 16:21]
               xx[1];					# 7
               xx[2] / (1 + xx[1]);		# 8
               xx[3];					# 9
               xx[4] * 10 * xx[7] / 10 / 9.5 ^ (-xx[6] - 1);	# 10
               xx[5] * 10 * xx[7] / 10 / 8.0 ^ (-xx[6] - 1);	# 11
               xx[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1];			# 12:17
               0;						# 18
               0;						# 19
               xx[12:13];				# 20:21
               xx[14]]					# 22

        b_boot[i, 1] = est[1]  # γ_s_shape
        b_boot[i, 2] = est[2]  # γ_s_on_09
        b_boot[i, 3] = est[3]  # γ_s_on_12
        b_boot[i, 4] = est[7]  # α
        b_boot[i, 5] = est[8]  # Δ_p_out
        b_boot[i, 6] = est[13] # (FIXED: r)
        b_boot[i, 7] = est[15] # c
        b_boot[i, 8] = est[4]  # σ_δ
        b_boot[i, 9] = est[12] # η
        b_boot[i,10] = est[9]                                      # (FIXED: γ_ns_shape)
        b_boot[i,11] =  est[10] .* est[9] # γ_ns_on_09
        b_boot[i,12] = (est[10] .* est[9] ./ 10) .^ (est[12] .+ 1) # (demand at p=10 2009 online)
        b_boot[i,13] =  est[11] .* est[9] # γ_ns_on_12
        b_boot[i,14] = (est[11] .* est[9] ./ 10) .^ (est[12] .+ 1) # (demand at p=10 2012 online)
        b_boot[i,15] = est[5]             # γ_ns_of_09_std
        b_boot[i,16] = (est[5] ./ 10) .^ (est[12] .+ 1)            # (demand at p=10 2009 offline)
        b_boot[i,17] = est[6] .*  est[5]  # γ_ns_of_09_loc (I CHANGED THIS to βlocal * gamma, used to be βlocal)
        b_boot[i,18] = est[16]            # γ_s_pop
        b_boot[i,19] = est[17]            # γ_ns_pop
        b_boot[i,20] = est[18] # (NOT USED: betacond)
        b_boot[i,21] = est[19] # (NOT USED: betapop)
        b_boot[i,22] = est[14]            # R_p
        b_boot[i,23] = est[20]            # s_R
        b_boot[i,24] = est[20] .* est[21] # μ_R
        b_boot[i,25] = 1.0 - est[22]      # R_q
    end

	# TODO: These statistics are computed but not returned anywhere
    betasigma_std_boot = [     std(b_boot[:,j])        for j=1:25]
    betasigma_boot_25  = [quantile(b_boot[:,j], 0.025) for j=1:25]
    betasigma_boot_5   = [quantile(b_boot[:,j], 0.05)  for j=1:25]
    betasigma_boot_10  = [quantile(b_boot[:,j], 0.1)   for j=1:25]
    betasigma_boot_90  = [quantile(b_boot[:,j], 0.9)   for j=1:25]
    betasigma_boot_95  = [quantile(b_boot[:,j], 0.95)  for j=1:25]
    betasigma_boot_975 = [quantile(b_boot[:,j], 0.975) for j=1:25]

    # Bootstrap welfare statistics
    boot_welfare = CSV.read(boot_out, DataFrame, header=false)
    boot_welfare = unique(boot_welfare)
    boot_welfare = boot_welfare[:, (end-8):end]

    welfare_std_boot = [     std(boot_welfare[:,j])        for j=1:9]
    welfare_boot_25  = [quantile(boot_welfare[:,j], 0.025) for j=1:9]
    welfare_boot_5   = [quantile(boot_welfare[:,j], 0.05)  for j=1:9]
    welfare_boot_10  = [quantile(boot_welfare[:,j], 0.1)   for j=1:9]
    welfare_boot_90  = [quantile(boot_welfare[:,j], 0.9)   for j=1:9]
    welfare_boot_95  = [quantile(boot_welfare[:,j], 0.95)  for j=1:9]
    welfare_boot_975 = [quantile(boot_welfare[:,j], 0.975) for j=1:9]

    # Note: column definitions are same as rows in Summary201609.xlsx
    if write_out
        CSV.write("$OUTPUT/bootstrap_estimates_$(vint).csv",         Tables.table(b_boot))
        CSV.write("$OUTPUT/bootstrap_welfare_estimates_$(vint).csv", boot_welfare)
    end
    return b_boot, boot_welfare
end

"""
Mapping of parameter vector to table names + locations:
. 1 => gamma0 shape parameter      == γ_s_shape   => 11
. 2 => E[gamma_0] 2009 online      == γ_s_on_09   => 8
. 3 => E[gamma_0] 2012 online      == γ_s_on_12   => 9
. 4 => alpha-1                     == α           => 12
. 5 => beta                        == Δ_p_out     => 13
X 6 => r                          (== r)
. 7 => c                           == c           => 14
. 8 => sigma_delta                 == σ_δ         => 7
. 9 => eta-1                       == η           => 6
X 10 => gamma_i shape parameter   (== γ_ns_shape)
. 11 => E[gamma_i] 2009 online     == γ_ns_on_09  => 2
  12 => demand at p=10 2009 online
. 13 => E[gamma_i] 2012 online     == γ_ns_on_12  => 3
  14 => demand at p=10 2012 online
. 15 => E[gamma_i] 2009 offline    == γ_ns_of_09_std
  16 => demand at p=10 2009 offline
. 17 => betalocal(*γ_ns_of_09_std) == γ_ns_of_09_loc
. 18 => lamda1                     == γ_s_pop     => 10
. 19 => lamda2                     == γ_ns_pop    => 5
X 20 => betacond
X 21 => betapop
. 22 => Pr(random price)           == R_p         => 15
. 23 => shape param. of rand price == s_R         => 17
. 24 => mean of random price       == μ_R         => 16
. 25 => Natural Disappear          == R_q         => 18

Welfare 2009 profit
Welfare 2009 CS nonshopper
Welfare 2009 CS shopper
Welfare 2012 profit
Welfare 2012 CS shopper
Welfare 2012 CS nonshopper
Welfare 2009 offline  profit
Welfare 2009 offline CS nonshopper
Welfare 2009 offline CS shopper
"""
function make_table_results(b_boot::Matrix{Float64}; table_title = "estimates.tex",
        table_rows = Dict([i => e for (i,e) in enumerate(
            [#=1=# (:γ_ns_of_09_std, "Mean arrival 2009 offline standard title (\$\\Bar{\\gamma}^{ns}_{09,of,std}\$)", "0.65 (0.17)"),
             #=2=# (:γ_ns_of_09_loc, "Mean arrival 2009 offline local interest (\$\\Bar{\\gamma}^{ns}_{12,of,std}\$)", "1.25 (0.37)"),
             #=3=# (:γ_ns_on_09, "Mean arrival 2009 online (\$\\Bar{\\gamma}^{ns}_{09,on}\$)", "14.90 (3.61)"),
             #=4=# (:γ_ns_on_12, "Mean arrival 2012 online (\$\\Bar{\\gamma}^{ns}_{12,on}\$)", "7.95 (2.00)"),
             #=5=# (:γ_ns_pop, "Popularity effect on arrival (\$\\gamma^{ns}_{Pop}\$)", "-1.36 (0.19)"),
             #=6=# (:η, "Nonshopper price elasticity (\$\\eta\$)", "1.87 (0.08)"),
             #"\\tp{Nonshopper condition elasticity (\$\\eta^c\$)}"),
             #=7=# (:σ_δ, "StDev of title-level unobservable (\$\\sigma_{\\delta}\$)", "1.16 (0.05)"),
             #=8=# (:γ_s_on_09, "Mean arrival 2009 online (\$\\Bar{\\gamma}^s_{09}\$)","5.65 (1.41)"),
             #=9=# (:γ_s_on_12, "Mean arrival 2012 online (\$\\Bar{\\gamma}^s_{12}\$)","14.86 (3.10)"),
             #=10=# (:γ_s_pop, "Popularity effect on arrival (\$\\gamma^{s}_{Pop}\$)","0.80 (0.09)"),
             #=11=# (:γ_s_shape, "Arrival distribution shape parameter (\$\\gamma^s_{shape}\$)","0.32 (0.03)"),
             #=12=# (:α, "Shopper price coefficient (\$\\alpha\$)", "15.77 (1.28)"),
             #=13=# (:Δ_p_out, "Outside option relative price (\$\\Delta p^{out}\$)", "0.16 (0.02)"),
             #=14=# (:c, "Effect of unobservable on price coefficient (\$c\$)", "-0.91 (0.07)"),
             #=15=# (:R_p, "Probability of randomly chosen price (\$R^p\$)", "0.26 (0.02)"),
             #=16=# (:μ_R, "Mean of randomly chosen price (\$\\mu^R\$)", "15.25 (0.80)"),
             #=17=# (:s_R, "Random price shape parameter (\$s^R\$)", "1.73 (0.09)"),
             #=18=# (:R_q, "Background disappearance rate (\$R^q\$)", "0.07 (0.01)")])]),
        label = "tab:my_label", include_paper = true)

    boot_mean = mean(b_boot, dims=1)
    boot_std  =  std(b_boot, dims=1)

    boot_ind = Dict([:γ_s_shape => 1, :γ_s_on_09 => 2, :γ_s_on_12 => 3, :α => 4,
					 :Δ_p_out => 5, # :r => 6,
                     :c => 7, :σ_δ => 8, :η => 9, # :γ_s_shape => 10,
                     :γ_ns_on_09 => 11, :γ_ns_on_12 => 13,
                     :γ_ns_of_09_std => 15, :γ_ns_of_09_loc => 17, #12, 14, 16 are demand at p=10
                     :γ_s_pop => 18, :γ_ns_pop => 19, # :βcond => 20, :βpop => 21
                     :R_p => 22, :s_R => 23, :μ_R => 24, :R_q => 25])

    io = open("tables/" * table_title, "w")
    if include_paper
        write(io, "\\begin{table}[h]\n\\centering\n\\begin{tabular}{c|c|c}")
        write(io, "\\toprule Parameter & New Est. & Paper Est. (SE)\\\\[2mm]")
    else
        write(io, "\\begin{table}[h]\n\\centering\n\\begin{tabular}{c|c}")
        write(io, "\\toprule Parameter & New Est. \\\\[2mm]")
    end

    for i=1:length(table_rows)
        θ_i = table_rows[i][1]
        if i==1
            write(io, "\\hline \n Nonshopper arrival and demand  \\\\ \n \\hline\n")
        elseif θ_i == :γ_s_on_09
            write(io, "\\hline \n Shopper arrival and demand  \\\\ \n \\hline\n")
        elseif θ_i == :R_p
            write(io, "\\hline \n Departures from fully rational model  \\\\ \n \\hline\n")
        end
        write(io, table_rows[i][2] * " & ")

        k = boot_ind[θ_i]
        boot_μ =  (θ_i == :η) ? 1 + boot_mean[k] : boot_mean[k]
        @printf(io, " %0.2f (%0.2f)", boot_μ, boot_std[k])

        include_paper ? @printf(io, " & %s", table_rows[i][3]) : nothing
        write(io, "\\\\ \n ")
    end
    write(io, "\\end{tabular}\\label{$label}\\end{table}")
    close(io)
end

import Base.zeros
zeros(x::Float64) = zeros(Int(round(x)))
zeros(x::Float64,y::Int64) = zeros(Int(round(x)), y)
zeros(x::Float64,y::Float64) = zeros(Int(round(x)), Int(round(y)))

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
  m, n = length(v1), length(v2)
  v1 = reshape(v1, m, 1)
  v2 = reshape(v2, 1, n)
  return repeat(v1, 1, n), repeat(v2, m, 1)
end

"""
```
function solveγPar(Dm, D0, dDm, dD0, γ0, p, r)
```
Works.
"""
function solveγPar(Dm, D0, dDm, dD0, γ0, p, r)
    A  = Dm .^ 2.0
    B  = (dDm .* r .* p) .+ (r .* Dm) .+ (2.0 .* Dm .* γ0 .* D0)
    C  = r .* p .* γ0 .* dD0 .+ r .* γ0 .* D0 .+ γ0 .* γ0 .* D0 .* D0
    γ1 = ((-B .+ ((B .^ 2 .- 4 .* A .* C) .+ 0im) .^ (0.5)) ./ (2 .* A))
    l1 = (real.(γ1) .> 0) .& (imag.(γ1) .≈ 0)
    return γ1, l1
end

function objective(x, x0, distpara0, γ0vec, δvec, data12, data09, bp)
    if x[3]<0 || x[4]<0 || x[12]>1 || x[12]< 0
        return Inf
    end
    xx = vcat(x[1:2], x0[3], x[3:5], x0[7], x[6:(length(x0)-2)])
    # if xx[4]<0 || xx[5]<0 || xx[14]>1 || xx[14]< 0
    #     return Inf
    # end
    return full_model(xx, distpara0, γ0vec, δvec, data12, data09, bp)
end

"""
```
function demandshopper(α, β, p, cdid, obs_w)
```
"""
function demandshopper(α, β, p, cdid, obs_w)
    expU = vec(exp.(p .* -α))
    sum1 = sum(sparse(collect(1:length(cdid)), vec(cdid), vec(obs_w) .* expU))'
            + exp(β .* α)

    f1 = expU ./ sum1
    f2 = -α .* expU .* (sum1 .- expU) ./ (sum1 .^ 2)
    f3 = α .^ 2 .* expU .* (sum1 .- expU) .* (sum1 .- 2 .* expU) ./ (sum1 .^ 3)

    replace!.([f1, f2, f3], NaN => 0)

    return f1, f2, f3, fill(sum1, size(cdid)), expU
end

"""
```
function obscalnewtest2015(βσ3, data, basellh, demandcal, p0, rounderr, WFcal)
```
Corresponds to Masao/obscalnewtest2015.m. Works!
"""
function obscalnewtest2015(βσ3, data, basellh, demandcal, p0::Vector{T},
                           rounderr; WFcal = false) where T<:Float64

    # [muγ0 α-1 β γishape γimean η-1 r lamda1 lamda2 βcond βpop βlocal olp δ c]
    numlist  = vec(data["numlist"])
    localint = vec(data["localint"])
    N        = Int(data["N"])
    M        = Int(data["M"])
    cdindex  = Int.(vec(data["cdindex"]))
    d_first  = Int.(vec(data["first"]))
    cond_dif = vec(data["conditiondif"])
    cdid     = vec(data["cdid"])
    obs_w    = vec(data["obsweight"])
    p        = vec(data["p"])

    γ0       = βσ3[1] .* (numlist .^ βσ3[8] ./ mean(numlist .^ βσ3[8]))
    α        = (βσ3[2] + 1) * βσ3[14] ^ βσ3[15]
    β        = βσ3[3] ./ βσ3[14] ^ βσ3[15]
    m        = βσ3[4]
    η        = βσ3[6] + 1
    r        = βσ3[7]
    βcond    = βσ3[10]
    olp      = βσ3[13]
    δ        = βσ3[14]
    γscale   = βσ3[5] ./ m .* (numlist .^ βσ3[9] ./ mean(numlist .^ βσ3[9])) .*
                   exp.(βσ3[12] .* localint)
    nat_disap = βσ3[16]

    ##########
    # Solve for demand and its 1st & 2nd order derivatives
    ##########
    D0, dD0, d2D0, sumpind, expU = demandshopper(α, β, p0 .- βcond .* cond_dif ./ α, cdid, obs_w)

    ## calculate γ1
    p1   = p .- rounderr
    Dm   = δ .* (p1 .^ (-η))
    dDm  = δ .* (-η) .* (p1 .^ (-η-1))
    d2Dm = δ .* (η) .* (η+1) .* (p1 .^ (-η-2))

    # Solve for the γ rationalizing price choice
    γ1, l1 = solveγPar(Dm, D0 .- rounderr .* dD0 .+ 0.5 .* d2D0 .* rounderr^2,
                       dDm, dD0 .- rounderr .* d2D0, γ0, p1, r)

    ## calculate γ2
    p2   = p .+ rounderr
    Dm   = δ .* (p2 .^ (-η))
    dDm  = δ .* (-η) .* (p2 .^ (-η-1))
    d2Dm = δ .* (η) .* (η+1) .* (p2 .^ (-η-2))

    # Solve for the γ rationalizing price choice
    γ2, l2 = solveγPar(Dm, D0 .+ rounderr .* dD0 .+ 0.5 .* d2D0 .* rounderr^2,
                       dDm, dD0 .+ rounderr .* d2D0, γ0, p2, r)

    SOC = r .* p2 .* (γ2 .* d2Dm .+ γ0 .* d2D0) .+ 2 .* (r .+ γ2 .* Dm .+ γ0 .*
            (D0 .+ rounderr .* dD0 .+ 0.5 .* d2D0 .* rounderr ^ 2)) .*
            (γ2 .* dDm .+ γ0 .* (dD0 .+ rounderr .* d2D0))

    γ1 = real.(γ1)
    γ2[imag.(γ2) .!= 0] .= 0
    γ2[real.(SOC) .> 0] .= 0
    γ2[real.(γ2)  .< 0] .= 0
    γ2 = real.(γ2)

    profit2 = p2 ./ (r ./ (γ2 .* Dm + γ0 .* D0) .+ 1)
    profitH = ((r .* (η-1) ./ δ) .^ (-1/η) ./ (1/(η-1)+1)) .* γ2 .^ (1/η)

    γ3, _ = solveγPar(Dm, 0, dDm, 0, 0, p2, r)
    γ2[real.(profitH) .> real.(profit2)] .= γ3[real.(profitH) .> real.(profit2)]

    γ1[real.(γ1) .< 0] .= .0
    γ1[real.(γ1) .> real.(γ2)] .= γ2[real.(γ1) .> real.(γ2)]

    γ1 = real.(γ1)
    γ2 = real.(γ2)

    Dm = δ .* (p .^ (-η))

    if demandcal == 1
        disap = data["disappear"]
        # D0 =demandshopper(α,β,p0- βcond.*data.cond_dif./α,data["cdid"],data["obsweight"]) #
        demandlh   = (disap > 0) - (2 .* disap - 1) .* exp(-0.166666667 .* (γ0.*D0 + (γ2 + γ1)./2 .* Dm)) .* nat_disap
        demandlhol = (disap > 0) - (2 .* disap - 1) .*  # Next line starts probability of nondisappear due to shopper demand
                        exp.(-0.166666667.*(γ0.*D0)) .* # Next line is due to nonshopper demand (taken expectation wrt to γi)
                        (1 .+ γscale .* 0.166666667 .* Dm) .^ -m .* nat_disap

        lip_o  = min((gamcdf(γ2, m, γscale) - gamcdf(γ1, m, γscale)), 1) .* demandlh .^ 3
        lip_ol = basellh .* demandlhol .^ 3 # Price likelihood. Next line starts disappear likelihood
    else
        γ_dist = Gamma.(m, γscale)
        lip_o  = min.([cdf(γ_dist[i], γ2[i]) - cdf(γ_dist[i], γ1[i]) for i=1:length(γ_dist)], 1) # TODO
        lip_ol = basellh # Price likelihood. Next line starts disappear likelihood.
    end

    liptemp = (1 - olp) .* lip_o + olp .* lip_ol
    olppost = vec(olp .* lip_ol ./ liptemp)
    lip = log.(liptemp)

    pi_v, CSns, CSs = zeros(N), zeros(N), zeros(N)

    if WFcal
        pi_v, CSns, CSs = welfaresimple(γ1, γ2, γscale .* m, γ0, olppost, Dm, D0, vec(p0),
                                        p, N, M, cdindex, d_first, vcat(α, β, η, r))
    end
    return lip, γ2, γ1, γ0, D0, Dm, pi_v, CSns, CSs
end

"""
```
function welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, pdif, data, scalarparas)
```
Corresponds to /Masao/welfaresimple.m. Works!
"""
function welfaresimple(γ1::Vector{T}, γ2::Vector{T}, γscale::Vector{T}, γ0::Vector{T},
                       olppost::Vector{T}, Dm::Vector{T}, D0::Vector{T}, pdif::Vector{T},
                       p::Vector{T}, N::S, M::S, cdindex::Vector{S}, d_first::Vector{S},
                       scalarparas::Vector{T}) where {S<:Int64, T<:Float64}
    α = scalarparas[1]
    β = scalarparas[2]
    η = scalarparas[3]
    r = scalarparas[4]

    γ1ave = 0.5 * (γ1 + γ2)

    pi_o  = p .* (γ1ave  .* Dm .+ γ0 .* D0) ./ (r .+ γ1ave  .* Dm .+ γ0 .* D0)
    pi_ol = p .* (γscale .* Dm .+ γ0 .* D0) ./ (r .+ γscale .* Dm .+ γ0 .* D0)
    pi_v  = olppost .* pi_ol .+ (1 .- olppost) .* pi_o

    CSns_o  = (1/(η-1)) .* γ1ave  .* Dm .* p ./ (r .+ γ1ave  .* Dm + γ0 .* D0)
    CSns_ol = (1/(η-1)) .* γscale .* Dm .* p ./ (r .+ γscale .* Dm + γ0 .* D0)
    CSns    = olppost .* CSns_ol .+ (1 .- olppost) .* CSns_o

    CSgain = zeros(N, 1)
    CSgain_OG = zeros(N, 1)

    N = 10000 # TODO: this looks like trouble

    mktsize = cdindex .- d_first .+ 1

    for k = 1:M
        ind_k    = d_first[k]:cdindex[k]

        best_p   = Vector{Float64}(undef, N)
        best_ind = Vector{Float64}(undef, N)
        β_p      = Vector{Float64}(undef, N)

        gumb_draw  = rand(rng, Gumbel(0,1), mktsize[k] + 1, N)

        rand_price = repeat(-[pdif[ind_k, 1] ; -β], 1, N) .- gumb_draw ./ α

        best, bestindex = vec.(findmax(rand_price, dims = 1))

        temp = sparse(vcat((x->x[2]).(bestindex), mktsize[k] + 1), 1:N+1, vcat(best .- rand_price[end,:], 1))

        CSgain_OG[d_first[k]:cdindex[k], 1] = sum(temp[1:mktsize[k],1:N], dims=2) ./ (sum(temp[1:mktsize[k],1:N] .> 0, dims=2) .+ 1e-5)
    end

    CSs_o  = γ0 .* D0 ./ (r .+ γ1ave  .* Dm .+ γ0 .* D0) .* CSgain
    CSs_ol = γ0 .* D0 ./ (r .+ γscale .* Dm .+ γ0 .* D0) .* CSgain
    CSs    = olppost .* CSs_ol + (1 .- olppost) .* CSs_o

    return pi_v, CSns, CSs
end

import Base.zeros
zeros(x::Float64) = zeros(Int(round(x)))
zeros(x::Float64,y::Int64) = zeros(Int(round(x)), y)
zeros(x::Float64,y::Float64) = zeros(Int(round(x)), Int(round(y)))

function ndgrid(v1::AbstractVector{T}, v2::AbstractVector{T}) where {T}
  m, n = length(v1), length(v2)
  v1 = reshape(v1, m, 1)
  v2 = reshape(v2, 1, n)
  (repeat(v1, 1, n), repeat(v2, m, 1))
end

function solveγPar(Dm, D0, dDm, dD0, γ0, p, r)

    A = Dm.*(Dm)
    B = (dDm) .* (r) .* (p) + r .* Dm + 2*Dm .* γ0 .* D0
    C = r .* p .* γ0 .* dD0 + r .* γ0 .* D0 + γ0 .* γ0 .* D0 .* D0

    γ1 = ((-B + (B .^ 2 .- 4 .* A .* C) .^ (0.5)) ./ (2 .* A))
    l1 = (γ1.re > 0) & (γ1.im == 0)

    return γ1, l1
end

function objective(x, x0, distpara0, γ0vec, deltavec, data12, data09, bp)
    if x[3]<0 || x[4]<0 || x[12]>1 || x[12]< 0
        return Inf
    end
    xx = vcat(x[1:2], x0[3], x[3:5], x0[7], x[6:(length(x0)-2)])
    # if xx[4]<0 || xx[5]<0 || xx[14]>1 || xx[14]< 0
    #     return Inf
    # end
    return full_model(xx, distpara0, γ0vec, deltavec, data12, data09, bp)
end

"""
```
function demandshopper(α, β, p, cdid, obsweight)
```
Runs!
"""
function demandshopper(α, β, p, cdid, obsweight)
    expU = vec(exp.(p .* -α))

    # sum utility for each book
    # sum1 = zeros(M)
    # sum1[1] = sum(expU(1:cdindex[1]))+exp(β)
    # for j=2:M
    #     sum1(j)=sum(expU(cdindex(j-1)+1:cdindex(j)))+exp(β)
    # end
    sum1 = sum(sparse(collect(1:length(cdid)), vec(cdid), vec(obsweight) .* expU))' + exp(β .* α[1])

    f1 = expU ./ sum1
    f2 = -α .* expU .* (sum1 .- expU) ./ (sum1 .^ 2)
    f3 = α .^ 2 .* expU .* (sum1 .- expU) .* (sum1 .- 2 .* expU) ./ (sum1 .^ 3)

    replace!.([f1, f2, f3], NaN => 0)

    return f1, f2, f3, fill(sum1, size(cdid)), expU
end

"""
```
function obscalnewtest2015(βsigma3, data, basellh, demandcal, p0, rounderr, WFcal)
```
Corresponds to Masao/obscalnewtest2015.m
"""
function obscalnewtest2015(βsigma3, data, basellh, demandcal, p0, rounderr, WFcal)

    # [muγ0 α-1 β γishape γimean η-1 r lamda1 lamda2 βcond βpop βlocal olp delta c]
    numlist      = data["numlist"]
    localint     = data["localint"]
    N            = Int(data["N"])
    conditiondif = data["conditiondif"]
    cdid         = vec(data["cdid"])
    obsweight    = vec(data["obsweight"])
    p            = data["p"]

    γ0 = fill(βsigma3[1], N, 1)' .* (numlist .^ βsigma3[8] ./ mean(numlist.^βsigma3[8]))
    α = fill((βsigma3[2]+1) * βsigma3[14] ^ βsigma3[15], N, 1)
    β = βsigma3[3] ./ βsigma3[14] ^ βsigma3[15]
    m = βsigma3[4]
    γscale = βsigma3[5] ./ m .* (numlist .^ βsigma3[9] ./ mean(numlist .^ βsigma3[9])) .* exp.(βsigma3[12].*localint)
    η = βsigma3[6]+1 # βsigma3(6*ones(N,1))'+1+βsigma3[11].*popular
    r = βsigma3[7]
    βcond = βsigma3[10]
    olp = βsigma3[13]
    delta = βsigma3[14]
    naturaldisappear = βsigma3[16]

    ##########
    # Solve for demand and its 1st & 2nd order derivatives
    ##########
    D0, dD0, d2D0, sumpind, expU = demandshopper(α, β, p0 - βcond .* conditiondif ./ α, cdid, obsweight)

    ## calculate γ1
    p = data["p"] - rounderr
    Dm=delta .* (p.^(-η))
    dDm=delta.*(-η).*(p.^(-η-1))
    d2Dm=delta.*(η).*(η+1).*(p.^(-η-2))

    #Solve for the γ that rationalize the price choice
    γ1, l1 = solveγPar(Dm,D0-rounderr.*dD0+0.5.*d2D0.*rounderr^2,dDm,dD0-rounderr.*d2D0,γ0,p,r)


    ## calculate γ2

    p = data["p"] + rounderr
    Dm=delta.*(p.^(-η))
    dDm=delta.*(-η).*(p.^(-η-1))
    d2Dm=delta.*(η).*(η+1).*(p.^(-η-2))

    #Solve for the γ that rationalize the price choice
    γ2, l2 = solveγPar(Dm,D0+rounderr.*dD0+0.5.*d2D0.*rounderr^2,dDm,dD0+rounderr.*d2D0,γ0,p,r)

    SOC = r.*p.*(γ2.*d2Dm + γ0.*d2D0) + 2*(r+γ2.*Dm + γ0.*(D0+rounderr.*dD0+0.5.*d2D0.*rounderr^2)).*(γ2.*dDm + γ0.*(dD0+rounderr.*d2D0))


    γ1[(im(γ1)~=0),1] = re(γ1((im(γ1)~=0),1))
    γ2[(im(γ2)~=0),1] = 0
    γ2[SOC .> 0, 1]  .= 0
    γ2 = maximum.(γ2, 0)

    profit2 = p./(r./(γ2.*Dm+γ0.*D0)+1)
    profitH = ((r * (η-1) / delta) ^ (-1/η) / (1/(η-1)+1)) .* γ2 .^ (1/η)
    γ3 = solveγPar(Dm,0,dDm,0,0,p,r)
    γ2[profitH .> profit2, 1] = γ3[profitH .> profit2,1]

    γ1 = min(max(γ1,0),γ2)

    Dm=delta.*(data.p.^(-η))

    ##
    if demandcal == 1

        # D0 =demandshopper(α,β,p0- βcond.*data.conditiondif./α,data["cdid"],data["obsweight"]) #
        demandlh = (data.disappear>0)-(2.*data.disappear -1).*exp(-0.166666667.*(γ0.*D0 + (γ2+γ1)./2.*Dm)).*naturaldisappear
        demandlhol = (data.disappear>0)-(2.*data.disappear - 1 ) ... Next line starts probability of nondisappear
            .*exp(-0.166666667.*(γ0.*D0)) ...due to shopper demand. Next line is due to nonshopper demand (taken expectation wrt to γi)
            .*(1 + γscale * 0.166666667.* Dm).^-m.*naturaldisappear

    # lip_o = interp1(pdfdata1,pdfdata2,γ1./γscale)./γscale.*dγidg.*demandlh.^3 #
    # lip_o = exp(-γ1./γscale)./γscale.*dγidg.*demandlh.^3 #
     lip_o = min((gamcdf(γ2,m,γscale) - gamcdf(γ1,m,γscale)),1).*demandlh.^3 #
        lip_ol = basellh.*  demandlhol.^3 #  # price likelihood. Next line starts disappear likelihood
    else
    # lip_o = interp1(pdfdata1,pdfdata2,γ1./γscale)./γscale.*dγidg #
    # lip_o = exp(-γ1./γscale)./γscale.*dγidg #
    lip_o = min(gamcdf(γ2,m,γscale) - gamcdf(γ1,m,γscale),1)
        lip_ol = basellh #  # price likelihood. Next line starts disappear likelihood
    end

    # lip_o(γ1<0,1) = 0
    # lip_o(lip_o<0,1) = 0
    liptemp = (1-olp).*lip_o + olp.*lip_ol
    olppost = olp.*lip_ol./liptemp
    lip = log(liptemp)

    pi_v = zeros(N)
    CSns = zeros(N)
    CSs  = zeros(N)

    if WFcal
        pi_v, CSns, CSs = welfaresimple(γ1, γ2, γscale.*m, γ0, olppost, Dm, D0, p0, data, [α[1] β η r])
    end
    return lip, γ2, γ1, γ0, D0, Dm, pi_v, CSns, CSs
end


"""
```
function welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, pdif, data, scalarparas)
```
Corresponds to /Masao/welfaresimple.m.
"""
function welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, pdif, data, scalarparas)

    α = scalarparas[1]
    β = scalarparas[2]
    η = scalarparas[3]
    r = scalarparas[4]

    γ1ave = 0.5 * (γ1 + γ2)

    pi_o  = data.p .* (γ1ave  .* Dm + γ0 .* D0) ./ (r + γ1ave  .* Dm + γ0 .* D0)
    pi_ol = data.p .* (γscale .* Dm + γ0 .* D0) ./ (r + γscale .* Dm + γ0 .* D0)
    pi_v  = olppost .* pi_ol + (1-olppost) .* pi_o

    CSns_o = (1/(η-1)) .* γ1ave .* Dm .* data.p ./ (r + γ1ave .* Dm + γ0 .* D0)
    CSns_ol = (1/(η-1)) .* γscale .* Dm .* data.p./(r + γscale .* Dm + γ0 .* D0)
    CSns = olppost .* CSns_ol + (1-olppost) .* CSns_o

    CSgain = zeros(data.N, 1)
    N = 10000

    cdindex = data["cdindex"]

    for k = 1:data.M

        mktsize = cdindex[k]-data["first"][k]+1
        randomprice = repeat(-[pdif(data["first"][k]:cdindex[k],1) ; -β],1,N) -
             evrnd(0,1,mktsize+1,N) ./ α[1]
        [best, bestindex] = max(randomprice)
        temp = sparse([bestindex mktsize+1],1:N+1,[best - randomprice(end,:) 1])
        CSgain[data["first"][k]:data["cdindex"][k],1] = sum(temp[1:mktsize,1:N],dims=2) ./ (sum(temp[1:mktsize,1:N] .> 0, dims=2) + 1e-5)

    end

    CSs_o  = γ0 .* D0 ./ (r + γ1ave  .* Dm + γ0 .* D0) .* CSgain
    CSs_ol = γ0 .* D0 ./ (r + γscale .* Dm + γ0 .* D0) .* CSgain
    CSs    = olppost .* CSs_ol + (1 - olppost) .* CSs_o
    return pi_v, CSns, CSs
end

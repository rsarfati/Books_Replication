# Based on file <fullmodelllhWFAug22newtest2015.m>
function full_model(x0, distpara0, γ0vec, δ_vec, data12, data09, bp; WFcal = false)
    # βσ5 = [γ0shape γ0mean09 γ0mean12 δ_σ
    #               γimeanbp βlocal alpha-1 β γishape  γimean09
    #               γimean12  eta-1 r olp c λ1 λ2 βcond βpop
    #               olm ol_θ naturaldisappear]
    # βσ4 = [γ0shape γ0mean09 γ0mean12 δ_σ γimeanbp
    #               alpha-1 β γishape  γimean09 γimean12  eta-1 r
    #               olp c λ1 λ2 βcond βpop βlocal olm ol_θ]
    βσ5 = [distpara0; x0[1]; x0[2]/(1+x0[1]); x0[3]; x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1);
           x0[5]*10*x0[7]/10/8^(-x0[6]-1); x0[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1]; 0; 0;
           x0[12:13]; x0[14]]

    rounderr = 0.025
    naturaldisappear = βσ5[22]

    βσ4 = βσ5[[1:5; 7:19; 6; 20; 21]]

    λ1    = βσ4[15]
    λ2    = βσ4[16]
    βcond = βσ4[17]
    βpop  = βσ4[18]
    olm   = βσ4[20]
    ol_θ  = βσ4[21]

    #TODO: ngrid very seldom necessary, needlessly memory intensive
    temp1, temp2 = ndgrid(γ0vec, δ_vec)
    γ0_δ_vec = hcat(vec(temp1), vec(temp2))

    Y    = length(temp1)
    M_09 = Int(data09["M"])
    N_09 = Int(data09["N"])
    M_12 = Int(data12["M"])
    N_12 = Int(data12["N"])

    first_09     = vecI64(data09["first"])
    p_09         = vecF64(data09["p"])
    pdif_09      = vecF64(data09["pdif"])
    numlist_09   = vecI64(data09["numlist"])
    cdindex_09   = vecI64(data09["cdindex"])
    cdid_09      = vecI64(data09["cdid"])

    cdid_12      = vecI64(data12["cdid"])
    first_12     = vecI64(data12["first"])
    p_12         = vecF64(data12["p"])
    pdif_12      = vecF64(data12["pdif"])
    numlist_12   = vecI64(data12["numlist"])
    cdindex_12   = vecI64(data12["cdindex"])
    obsweight_09 = vecF64(data09["obsweight"])
    obsweight_12 = vecF64(data12["obsweight"])

    basellh_09   = pdf.(Gamma(olm, ol_θ), p_09) * 2 * rounderr

    bp_p = vecF64(bp["p"])

    # Calculation for 09 data
    ltot_09 = zeros(M_09)
    llhβ_09 = zeros(M_09, Y)  # Record the log likelihood at a fixed β for a title

    # TODO: this gets run a LOT, present run time is...
    # 2.255119 seconds (59.73 k allocations: 1.512 GiB, 5.42% gc time) per call
    println("Iterating for γ0...")
    out_09 = if parallel
       @distributed (hcat) for i in 1:Y
           obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6, 7, 8, 9, 11, 12]]; λ1; λ2;
                              βcond; βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                             data09, basellh_09, pdif_09, rounderr;
                             demandcal = false, WFcal = WFcal)
       end
    else
       hcat([obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6, 7, 8, 9, 11, 12]]; λ1; λ2;
                                βcond; βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                               data09, basellh_09, pdif_09, rounderr;
                               demandcal = false, WFcal = WFcal) for i=1:Y]...)
    end
    lip_09  = out_09[0*N_09+1:1*N_09,:] # likelihood of each observation at each β
    γ2_09   = out_09[1*N_09+1:2*N_09,:]
    γ1_09   = out_09[2*N_09+1:3*N_09,:]
    γ0_09   = out_09[3*N_09+1:4*N_09,:]
    D0_09   = out_09[4*N_09+1:5*N_09,:]
    Dm_09   = out_09[5*N_09+1:6*N_09,:]
    pi_09   = out_09[6*N_09+1:7*N_09,:]
    CSns_09 = out_09[7*N_09+1:8*N_09,:]
    CSs_09  = out_09[8*N_09+1:9*N_09,:]

    for k = 1:M_09
        llhβ_09[k,:] .= sum(lip_09[first_09[k]:cdindex_09[k],:])
    end
    maxtemp_09 = maximum(llhβ_09, dims=2)
    llhadj_09  = exp.(llhβ_09 - repeat(maxtemp_09, 1, Y))

    ## Calculation for 2012 data
    ltot_12    = zeros(M_12)
    llhβ_12    = zeros(M_12, Y)  # record the log likelihood at a fixed β for a title
    basellh_12 = pdf.(Gamma(olm, ol_θ), p_12) * 2 * rounderr

    println("Iterating for γ0...")
    out_12 = if parallel
       @distributed (hcat) for i in 1:Y
           obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6, 7, 8, 10, 11, 12]]; λ1; λ2;
                              βcond; βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; naturaldisappear],
                             data12, basellh12, pdif_12, rounderr;
                             demandcal = true, WFcal = WFcal)
       end
    else
       hcat([obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6, 7, 8, 10, 11, 12]]; λ1; λ2;
                                βcond; βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; naturaldisappear],
                               data12, basellh12, pdif_12, rounderr;
                               demandcal = true, WFcal = WFcal) for i=1:Y]...)
    end
    lip_12  = out_12[0*N_12+1:1*N_12,:]
    γ2_12   = out_12[1*N_12+1:2*N_12,:]
    γ1_12   = out_12[2*N_12+1:3*N_12,:]
    γ0_12   = out_12[3*N_12+1:4*N_12,:]
    D0_12   = out_12[4*N_12+1:5*N_12,:]
    Dm_12   = out_12[5*N_12+1:6*N_12,:]
    pi_12   = out_12[6*N_12+1:7*N_12,:]
    CSns_12 = out_12[7*N_12+1:8*N_12,:]
    CSs_12  = out_12[8*N_12+1:9*N_12,:]

    for k=1:M_12
        llhβ_12[k,:] .= sum(lip_12[first_12[k]:cdindex_12[k],:])
    end

    maxtemp_12 = maximum(llhβ_12, dims = 2)
    llhadj_12  = exp.(llhβ_12 .- repeat(maxtemp_12,1,Y))

    ## Calculation for 09 offline data
    basellhb = pdf.(Gamma(olm, ol_θ), bp_p) * 2 * rounderr

    getbmean(γ_l) = -sum(obscalnewtest2015([0.; βσ4[[6, 7, 8]]; abs(γ_l[1]);
                                            βσ4[[11, 12]]; λ1; λ2; βcond;
                                            βpop; γ_l[2]; βσ4[13]; 1.; 1.; 1.],
                                            bp, basellhb, bp_p,
                                            rounderr; demandcal = false,
                                            WFcal = false)[1:length(basellhb),:])

    function integγ0(γinput; return_all = false)
        γ0shape = γinput[1]
        γ0_θ_09 = γinput[2] / γinput[1]
        γ0_θ_12 = γinput[3] / γinput[1]
        δ_mean  = -0.5 * γinput[4] ^ 2
        δ_σ     = γinput[4]

        # Calculate the importance of the grid points
        γ0_cdf_vec = vcat(-cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec[1]),
                         cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec),
                         2 .- cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec[end]))
        imp1 = (γ0_cdf_vec[3:end] - γ0_cdf_vec[1:end-2]) ./ 2

        γ0_cdf_vec = vcat(-cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec[1]),
                         cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec),
                         2 .- cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec[end]))

        imp2 = (γ0_cdf_vec[3:end] - γ0_cdf_vec[1:end-2]) ./ 2

        δ_cdf_vec = [    -cdf( Normal(δ_mean, δ_σ), log( δ_vec[1]));
                          cdf.(Normal(δ_mean, δ_σ), log.(δ_vec));
                     2 .- cdf.(Normal(δ_mean, δ_σ), log( δ_vec[end]))]
        imp3 = (δ_cdf_vec[3:end] - δ_cdf_vec[1:end-2]) ./ 2

        imp_09 = imp1 * imp3'
        imp_12 = imp2 * imp3'

        ltot_09 = maxtemp_09 + log.(llhadj_09 * vec(imp_09))
        ltot_12 = maxtemp_12 + log.(llhadj_12 * vec(imp_12))

        if return_all
            return imp_09, imp_12, ltot_09, ltot_12
        end
        return -(sum(ltot_09) + sum(ltot_12))
    end

    res = optimize(integγ0, βσ4[1:4])
    distpara1, f1 = res.minimizer, res.minimum

    res = optimize(getbmean, βσ5[5:6])
    distpara2, f2 = res.minimizer, res.minimum

    out_bp = obscalnewtest2015(vcat(0, βσ4[[6, 7, 8]],
        abs(distpara2[1]), βσ4[[11, 12]], λ1, λ2, βcond, βpop, distpara2[2],
        βσ4[13], 1, 1, 1), bp, basellhb, bp_p, rounderr; demandcal = false, WFcal = WFcal)

    N_bp  = length(bp_p)
    WF_bp = zeros(N_bp, 3)

    lipb    = out_bp[0*N_bp+1:1*N_bp,:]
    γ2_bp   = out_bp[1*N_bp+1:2*N_bp,:]
    γ1_bp   = out_bp[2*N_bp+1:3*N_bp,:]
    γ0_bp   = out_bp[3*N_bp+1:4*N_bp,:]
    D0_bp   = out_bp[4*N_bp+1:5*N_bp,:]
    Dm_bp   = out_bp[5*N_bp+1:6*N_bp,:]
    WF_bp[:, 1] .= vec(out_bp[6*N_bp+1:7*N_bp,:])
    WF_bp[:, 2] .= vec(out_bp[7*N_bp+1:8*N_bp,:])
    WF_bp[:, 3] .= vec(out_bp[8*N_bp+1:9*N_bp,:])

    f = f1 + f2
    distpara = [distpara1; distpara2]

    if WFcal
        imp_09, imp_12, ltot_09, ltot_12 = integγ0(distpara1; return_all = true)
        WF_09       = zeros(N_09, 3)
        WF_12       = zeros(N_12, 3)
        AveWF_09    = zeros(M_09, 3)
        AveWF_12    = zeros(M_12, 3)
        BestVals_09 = zeros(N_09,13)
        BestVals_12 = zeros(N_12,13)

        for k = 1:M_09
            RPpost = llhadj_09[k,:]' .* vec(imp_09) ./ exp.(ltot_09[k] - maxtemp_09[k])
            ind_k  = first_09[k]:cdindex_09[k]

            WF_09[ind_k,1]  = pi_09[ind_k,:]   * RPpost
            WF_09[ind_k,2]  = CSns_09[ind_k,:] * RPpost
            WF_09[ind_k,3]  = CSs_09[ind_k,:]  * RPpost
            obsweight        = obsweight_09[ind_k, 1] ./ sum(obsweight_09[ind_k,1])
            AveWF_09[k,1:3] = obsweight' * WF_09[ind_k,1:3]
            y_max            = argmax(llhadj_09[k,:])

            BestVals_09[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max, :], llhβ_09[k, y_max] ./ length(ind_k)), 1, length(ind_k))',
                exp.(lip_09[ind_k, y_max]), basellh_09[ind_k,1], γ0_09[ind_k,y_max], γ1_09[ind_k,y_max], γ2_09[ind_k,y_max],
                Dm_09[ind_k,y_max], D0_09[ind_k,y_max], pi_09[ind_k,y_max], CSns_09[ind_k,y_max], CSs_09[ind_k,y_max])

            RPpost = llhadj_12[k,:] .* vec(imp_12) ./ exp(ltot_12[k]-maxtemp_12[k])

            ind_k           = first_12[k]:cdindex_12[k]
            WF_12[ind_k,1]  = pi_12[ind_k,:]  * RPpost
            WF_12[ind_k,2]  = CSns_12[ind_k,:]* RPpost
            WF_12[ind_k,3]  = CSs_12[ind_k,:] * RPpost
            obsweight       = obsweight_12[ind_k,1] ./ sum(obsweight_12[ind_k,1])
            AveWF_12[k,1:3] = obsweight' * WF_12[ind_k, 1:3]

            y_max = argmax(llhadj_12[k,:])
            BestVals_12[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max,:], llhβ_12[k,y_max] ./ length(ind_k))', length(ind_k), 1),
                exp.(lip_12[ind_k, y_max]), basellh_12[ind_k,1],
                γ0_12[ind_k,y_max], γ1_12[ind_k,y_max], γ2_12[ind_k,y_max],
                Dm_12[ind_k,y_max], D0_12[ind_k,y_max], pi_12[ind_k,y_max], CSns_12[ind_k,y_max], CSs_12[ind_k,y_max])
        end

        fWF = Dict{String,Any}()
        fWF["BestVals_09"] = hcat(cdid_09, numlist_09, p_09, pdif_09, BestVals_09)
        fWF["BestVals_12"] = hcat(cdid_12, numlist_12, p_12, pdif_12, BestVals_12)
        fWF["BestValsbp"] = hcat(bp["cdid"], bp["numlist"], bp["p"],
                                 repeat([0.0, 0.0, 1.0, 1.0]', Int.(bp["N"]), 1),
                                 lipb, basellhb, γ0_bp, γ1_bp, γ2_bp, Dm_bp, D0_bp, WF_bp)
        fWF["AveWF_09"] = AveWF_09
        fWF["AveWF12"]  = AveWF_12
        fWF["WF_09"]    = WF_09
        fWF["WF12"]     = WF_12
        fWF["WFbp"]     = WF_bp
    end

    fother = Dict{String,Any}()
    fother["lip12"]    = lip_12
    fother["llhβ12"]   = llhβ_12
    fother["lip_09"]   = lip_09
    fother["llhβ_09"]  = llhβ_09
    fother["γ0_δ_vec"] = γ0_δ_vec
    fother["imp_09"]   = imp_09
    fother["imp_12"]   = imp_12
    fother["ltot_09"]  = ltot_09
    fother["ltot12"]   = ltot_12
    fother["imp_09"]   = imp_09

    fother["γ1_09"] = γ1_09
    fother["γ112"]  = γ1_12

    fother["lipb"] = lipb
    fother["γ1bp"] = γ1_bp

    return f, distpara, fother, fWF
end


function objective(x, x0, distpara0, γ0vec, δvec, data12, data09, bp)
    if x[3]<0 || x[4]<0 || x[12]>1 || x[12]< 0
        return Inf
    end
    xx = vcat(x[1:2], x0[3], x[3:5], x0[7], x[6:(length(x0)-2)])

    #TODO: run this at some point
    @assert x[3] == xx[4]
    @assert x[4] == xx[5]
    @assert x[12] == xx[14]

    # if xx[4]<0 || xx[5]<0 || xx[14]>1 || xx[14]< 0
    #     return Inf
    # end

    return full_model(xx, distpara0, γ0vec, δvec, data12, data09, bp)
end

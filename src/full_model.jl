# Based on file <fullmodelllhWFAug22newtest2015.m>
function full_model(x0, distpara0, γ0vec, δ_vec, data12, data09, bp; WFcal = false)
    # βσ5 = [γ0shape γ0mean09 γ0mean12 δ_σ
    #               γimeanbp βlocal alpha-1 β γishape  γimean09
    #               γimean12  eta-1 r olp c λ1 λ2 βcond βpop
    #               olm ol_θ naturaldisappear]
    # βσ4 = [γ0shape γ0mean09 γ0mean12 δ_σ γimeanbp
    #               alpha-1 β γishape  γimean09 γimean12  eta-1 r
    #               olp c λ1 λ2 βcond βpop βlocal olm ol_θ]
    βσ5 = vcat(distpara0, x0[1], x0[2]/(1+x0[1]), x0[3],
        x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1), x0[5]*10*x0[7]/10/8^(-x0[6]-1),
        x0[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1], 0, 0, x0[12:13], x0[14])

    rounderr = 0.025
    naturaldisappear = βσ5[22]

    βσ4 = βσ5[vcat(1:5, 7:19, 6, 20, 21)]

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
    first_09 = Int.(vec(data09["first"]))
    p_09 = vec(data09["p"])
    pdif_09 = vec(data09["pdif"])
    numlist_09 = Int.(vec(data09["numlist"]))
    cdindex_09 = Int.(vec(data09["cdindex"]))
    cdid_09 = Int.(vec(data09["cdid"]))

    cdid_12 = Int.(vec(data12["cdid"]))
    first_12 = Int.(vec(data12["first"]))
    p_12 = vec(data12["p"])
    pdif_12 = vec(data12["pdif"])
    numlist_12 = Int.(vec(data12["numlist"]))
    cdindex_12 = Int.(vec(data12["cdindex"]))
    obsweight_09 = vec(data09["obsweight"])
    obsweight_12 = vec(data12["obsweight"])

    ## Calculation for 09 data
    ltot_09 = zeros(M_09)
    llhβ_09 = zeros(M_09, Y)  # record the log likelihood at a fixed β for a title
    lip_09  = zeros(N_09, Y)  # likelihood of each observation at each β
    γ1_09, γ2_09, γ0_09, D0_09, Dm_09 = lip_09, lip_09, lip_09, lip_09, lip_09

    basellh_09 = pdf.(Gamma(olm, ol_θ), p_09) * 2 * rounderr

    pi_09, CSns_09, CSs_09 = zeros(N_09, Y), zeros(N_09, Y), zeros(N_09, Y)

    # Iterate for γ0
    # TODO: this gets run a LOT, present run time:
    # 2.255119 seconds (59.73 k allocations: 1.512 GiB, 5.42% gc time)
    #@parallel
    for i = 1:Y
        @show i
        lip_09[:,i], γ2_09[:,i], γ1_09[:,i], γ0_09[:,i], D0_09[:,i], Dm_09[:,i],
        pi_09[:,i], CSns_09[:,i], CSs_09[:,i] = obscalnewtest2015(
                     vcat(γ0_δ_vec[i, 1], βσ4[[6, 7, 8, 9, 11, 12]], λ1, λ2,
                          βcond, βpop, 0, βσ4[13], γ0_δ_vec[i, 2], βσ4[14], 1),
                     data09, basellh_09, pdif_09, rounderr; demandcal = false, WFcal = WFcal)
    end
    for k = 1:M_09
        llhβ_09[k,:] .= sum(lip_09[first_09[k]:cdindex_09[k],:])
    end
    maxtemp_09 = maximum(llhβ_09, dims=2)
    llhadj_09  = exp.(llhβ_09 - repeat(maxtemp_09, 1, Y))

    ## Calculation for 2012 data
    lip12  = zeros(N_12, Y)
    γ112   = lip12
    γ212   = lip12
    γ012   = lip12
    D012   = lip12
    Dm12   = lip12
    ltot12 = zeros(M_12)
    llhβ12 = zeros(M_12, Y)  # record the log likelihood at a fixed β for a title
    basellh12 = pdf.(Gamma(olm, ol_θ), p_12) * 2 * rounderr
    pi12, CSns12, CSs12 = zeros(N_12, Y), zeros(N_12, Y), zeros(N_12, Y)

    # iterate for βs
    #@parallel
    for i = 1:Y
        @show i
        lip12[:,i], γ212[:,i], γ112[:,i], γ012[:,i], D012[:,i],
        Dm12[:,i], pi12[:,i], CSns12[:,i], CSs12[:,i] = obscalnewtest2015(
            vcat(γ0_δ_vec[i,1], βσ4[[6, 7, 8, 10, 11, 12]], λ1, λ2,
            βcond, βpop, 0, βσ4[13], γ0_δ_vec[i, 2], βσ4[14],
            naturaldisappear), data12, basellh12, 1, pdif_12, rounderr; WFcal = WFcal)
    end
    for k=1:M_12
        llhβ12[k,:] .= sum(lip12[first_12[k]:cdindex_12[k],:])
    end

    maxtemp_12 = maximum(llhβ12,dims=2)
    llhadj_12  = exp.(llhβ12 .- repeat(maxtemp12,1,Y))

    ## Calculation for 09 offline data
    basellhb = pdf.(Gamma(olm, ol_θ), bp["p"]) * 2 * rounderr

    getbmean(γ_l) = -sum(obscalnewtest2015(vcat(0, βσ4[[6, 7, 8]], abs(γ_l[1]),
                                                βσ4[[11, 12]], λ1, λ2, βcond,
                                                βpop, γ_l[2], βσ4[13], 1, 1, 1),
                                                bp, basellhb, 0, vec(bp["p"]),
                                                rounderr, false)[1])

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

        # TODO: find every use of quantile and verify it wasn't meant to be CDF
        δ_cdf_vec = vcat(-cdf(Normal(δ_mean, δ_σ), log(δ_vec[1])),
                        cdf.(Normal(δ_mean, δ_σ), log.(δ_vec)),
                        2 .- cdf.(Normal(δ_mean,δ_σ), log(δ_vec[end])))
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

    WFbp = zeros(length(lipb), 3)
    lipb, γ2bp, γ1bp, γ0bp, D0bp, Dmbp,
    WFbp[:,1], WFbp[:,2], WFbp[:,3] = obscalnewtest2015(vcat(0, βσ4[[6, 7, 8,]],
        abs(distpara2[1]), βσ4[[11, 12]], λ1, λ2, βcond, βpop, distpara2[2],
        βσ4[13], 1, 1, 1), bp, basellhb, 0, bp["p"], rounderr; WFcal = WFcal)

    f = f1 + f2
    distpara = vcat(distpara1, distpara2)

    if WFcal
        imp_09, imp_12, ltot_09, ltot_12 = integγ0(distpara1; return_all = true)
        WF_09       = zeros(N_09, 3)
        WF12        = zeros(N_12, 3)
        AveWF_09    = zeros(M_09, 3)
        AveWF12     = zeros(M_12, 3)
        BestVals_09 = zeros(N_09,13)
        BestVals_12 = zeros(N_12,13)

        for k = 1:M_09
            RPpost = llhadj_09[k,:] .* vec(imp_09) ./ exp(ltot_09[k] - maxtemp_09[k])
            ind_k  = first_09[k]:cdindex_09[k]

            WF_09[ind_k, 1]  = pi_09[ind_k,:]   * RPpost
            WF_09[ind_k, 2]  = CSns_09[ind_k,:] * RPpost
            WF_09[ind_k, 3]  = CSs_09[ind_k,:]  * RPpost
            obsweight        = obsweight_09[ind_k, 1] ./ sum(obsweight_09[ind_k,1])
            AveWF_09[k, 1:3] = obsweight' * WF_09[ind_k, 1:3]
            y_max            = argmax(llhadj_09[k,:])

            BestVals_09[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max, :], llhβ_09[k, y_max] ./ length(ind_k)), 1, length(ind_k))',
                exp.(lip_09[ind_k, y_max]), basellh_09[ind_k,1], γ0_09[ind_k,y_max], γ1_09[ind_k,y_max], γ2_09[ind_k,y_max],
                Dm_09[ind_k,y_max], D0_09[ind_k,y_max], pi_09[ind_k,y_max], CSns_09[ind_k,y_max], CSs_09[ind_k,y_max])

            RPpost = llhadj12[k,:] .* vec(imp_12) ./ exp(ltot_12[k]-maxtemp12[k])

            ind_k          = first_12[k]:cdindex_12[k]
            WF12[ind_k,1]  = pi12[ind_k,:]  * RPpost
            WF12[ind_k,2]  = CSns12[ind_k,:]* RPpost
            WF12[ind_k,3]  = CSs12[ind_k,:] * RPpost
            obsweight      = obsweight_12[ind_k,1] ./ sum(obsweight_12[ind_k,1])
            AveWF12[k,1:3] = obsweight' * WF12[ind_k, 1:3]

            y_max = argmax(llhadj12[k,:])
            BestVals_12[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max,:], llhβ12[k,y_max] ./ length(ind_k))', length(ind_k), 1),
                exp.(lip12[ind_k, y_max]), basellh12[ind_k,1],
                γ012[ind_k,y_max], γ112[ind_k,y_max], γ212[ind_k,y_max],
                Dm12[ind_k,y_max], D012[ind_k,y_max], pi12[ind_k,y_max], CSns12[ind_k,y_max], CSs12[ind_k,y_max])
        end

        fWF = Dict{String,Any}()
        fWF["BestVals_09"] = hcat(cdid_09, numlist_09, p_09, pdif_09, BestVals_09)
        fWF["BestVals_12"] = hcat(cdid_12, numlist_12, p_12, pdif_12, BestVals_12)
        fWF["BestValsbp"] = hcat(bp["cdid"], bp["numlist"], bp["p"],
                                 repeat([0.0, 0.0, 1.0, 1.0]', Int.(bp["N"]), 1),
                                 lipb, basellhb, γ0bp, γ1bp, γ2bp, Dmbp, D0bp, WFbp)
        fWF["AveWF_09"] = AveWF_09
        fWF["AveWF12"]  = AveWF12
        fWF["WF_09"]    = WF_09
        fWF["WF12"]     = WF12
        fWF["WFbp"]     = WFbp
    end

    fother = Dict{String,Any}()
    fother["lip12"]    = lip12
    fother["llhβ12"]   = llhβ12
    fother["lip_09"]   = lip_09
    fother["llhβ_09"]  = llhβ_09
    fother["γ0_δ_vec"] = γ0_δ_vec
    fother["imp_09"]   = imp_09
    fother["imp_12"]   = imp_12
    fother["ltot_09"]  = ltot_09
    fother["ltot12"]   = ltot12
    fother["imp_09"]   = imp_09

    fother["γ1_09"] = γ1_09
    fother["γ112"]  = γ112

    fother["lipb"] = lipb
    fother["γ1bp"] = γ1bp

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

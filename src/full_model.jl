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
    M_09 = data09["M"]
    N_09 = data09["N"]
    M_12 = data12["M"]
    N_12 = data12["N"]

    ## Calculation for 09 data
    ltot_09 = zeros(M_09)
    llhβ_09 = zeros(M_09, Y)  # record the log likelihood at a fixed β for a title
    lip_09  = zeros(N_09, Y)  # likelihood of each observation at each β
    γ1_09, γ2_09, γ0_09, D0_09, Dm_09 = lip_09, lip_09, lip_09, lip_09, lip_09
    basellh_09 = pdf.(Gamma(olm, ol_θ), data09["p"]) * 2 * rounderr

    pi_09, CSns_09, CSs_09 = zeros(N_09, Y), zeros(N_09, Y), zeros(N_09, Y)

    # Iterate for γ0
    #@parallel
    for i = 1#:Y
        obscalnewtest2015(vcat(γ0_δ_vec[i, 1],
                     βσ4[[6, 7, 8, 9, 11, 12]], λ1, λ2, βcond, βpop, 0, βσ4[13], γ0_δ_vec[i, 2],
                     βσ4[14], 1), data09, basellh_09, 0, data09["pdif"], rounderr, WFcal)

        lip_09[:,i], γ2_09[:,i], γ1_09[:,i], γ0_09[:,i], D0_09[:,i], Dm_09[:,i],
        pi_09[:,i], CSns_09[:,i], CSs_09[:,i] = obscalnewtest2015(vcat(γ0_δ_vec[i, 1],
                     βσ4[[6, 7, 8, 9, 11, 12]], λ1, λ2, βcond, βpop, 0, βσ4[13], γ0_δ_vec[i, 2],
                     βσ4[14], 1), data09, basellh_09, 0, data09["pdif"], rounderr, WFcal)
    end
    lip_09, γ2_09, γ1_09, γ0_09, D0_09, Dm_09, pi_09, CSns_09, CSs_09 = obscalnewtest2015([γ0_δ_vec[:,1] βσ4[[6 7 8 9 11 12]] λ1 λ2 βcond βpop 0 βσ4[13] γ0_δ_vec[:,2] βσ4[14] 1],data09,basellh_09,0,data09["pdif"],rounderr,WFcal)

    for k=1:M_09
        llhβ_09[k, :] = sum(lip_09(data09["first"][k]:data09["cdindex"][k],:))
    end
    maxtemp_09 = max(llhβ_09,[],2)
    llhadj_09  = exp(llhβ_09 - repmat(maxtemp_09,1,Y))

    ## Calculation for 2012 data
    lip12     = zeros(N_12, Y)
    γ112  = lip12
    γ212  = lip12
    γ012  = lip12
    D012      = lip12
    Dm12      = lip12
    ltot12    = zeros(M_12)
    llhβ12 = zeros(M_12, Y)  # record the log likelihood at a fixed β for a title
    basellh12 = pdf(Gamma(olm, ol_θ), data12["p"]) * 2 * rounderr
    pi12, CSns12, CSs12 = zeros(N_12, Y), zeros(N_12, Y), zeros(N_12, Y)

    # iterate for βs
    @parallel for i = 1:Y
        lip12[:,i], γ212[:,i], γ112[:,i], γ012[:,i], D012[:,i],
        Dm12[:,i], pi12[:,i], CSns12[:,i], CSs12[:,i] = obscalnewtest2015(
            vcat(γ0_δ_vec[i,1], βσ4[[6 7 8 10 11 12]], λ1, λ2,
            βcond, βpop, 0, βσ4[13], γ0_δ_vec[i,2], βσ4[14],
            naturaldisappear), data12, basellh12, 1, data12["pdif"], rounderr, WFcal)
    end
    for k=1:M_12
        llhβ12[k,:] = sum(lip12(data12["first"][k]:data12["cdindex"][k],:))
    end

    maxtemp12 = max(llhβ12,[],2)
    llhadj12  = exp(llhβ12 - fill(maxtemp12,1,Y))

    ## Calculation for 09 offline data
    basellhb = pdf(Gamma(olm, ol_θ), bp["p"]) * 2 * rounderr

    getbmean(γ_l) = -sum(obscalnewtest2015(vcat(0, βσ4[[6 7 8]], abs(γ_l[1]),
                                                βσ4[[11, 12]], λ1, λ2, βcond,
                                                βpop, γ_l[2], βσ4[13], 1, 1, 1),
                                           bp, basellhb, 0, bp["p"], rounderr, 0))

    function integγ0(γinput)
        γ0shape   = γinput[1]
        γ0_θ_09 = γinput[2]/γinput[1]
        γ0_θ_12 = γinput[3]/γinput[1]
        δ_mean     = -0.5*γinput[4]^2
        δ_σ    = γinput[4]

        # Calculate the imp of the grid points
        γ0cdfvec = [-gamcdf(γ0vec[1],γ0shape,γ0_θ_09) gamcdf(γ0vec,γ0shape,γ0_θ_09) (2 - gamcdf(γ0vec[end],γ0shape,γ0_θ_09))]
        imp1 = (γ0cdfvec[3:end]-γ0cdfvec[1:end-2])/2
        γ0cdfvec = [-gamcdf(γ0vec[1],γ0shape,γ0_θ_12) gamcdf(γ0vec,γ0shape,γ0_θ_12) (2 - gamcdf(γ0vec[end],γ0shape,γ0_θ_12))]
        imp2 = (γ0cdfvec[3:end]-γ0cdfvec[1:end-2])/2
        δ_cdfvec = [-normcdf(log(δ_vec[1]),δ_mean,δ_σ) normcdf(log(δ_vec),δ_mean,δ_σ) (2 - normcdf(log(δ_vec[end]),δ_mean,δ_σ))]
        imp3 = (δ_cdfvec[3:end]-δ_cdfvec[1:end-2])/2
        imp_09 = imp1' * imp3
        imp12 = imp2' * imp3

        ltot_09 = maxtemp_09 + log.(llhadj_09 .* vec(imp_09))
        ltot12 = maxtemp12 + log.(llhadj12 .* vec(imp12))
        return -(sum(ltot_09) + sum(ltot12))
    end

    [distpara1,f1,~,fmindisplay] = fminunc(@integγ0, βσ4[1:4], optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    disp(fmindisplay.funcCount)

    [distpara2,f2] = fminunc(@getbmean,βσ5[5:6],optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    lipb,γ2bp,γ1bp,γ0bp,D0bp,Dmbp,WFbp[:,1],WFbp[:,2],WFbp[:,3] = obscalnewtest2015([0 βσ4([6 7 8]) abs(distpara2[1]) βσ4([ 11 12 ]) λ1 λ2 βcond βpop distpara2[2] βσ4[13] 1 1 1],bp,basellhb,0,bp["p"],rounderr,WFcal)

    f = f1 + f2
    distpara = [distpara1, distpara2]

    if WFcal
        WF_09       = zeros(N_09, 3)
        WF12       = zeros(N_12, 3)
        AveWF_09    = zeros(M_09, 3)
        AveWF12    = zeros(M_12, 3)
        BestVals_09 = zeros(N_09,13)
        BestVals12 = zeros(N_12,13)

        for k = 1:M_09
            RPpost = llhadj_09[k,:]' .* imp_09[:]/exp(ltot_09[k] - maxtemp_09[k])
            ind_k  = data09["first"][k]:data09["cdindex"][k]

            WF_09[ind_k, 1]  = pi_09[ind_k,:]*RPpost
            WF_09[ind_k, 2]  = CSns_09[ind_k,:]*RPpost
            WF_09[ind_k, 3]  = CSs_09[ind_k,:]*RPpost
            obsweight       = data09["obsweight"][ind_k, 1] ./ sum(data09["obsweight"][ind_k, 1])
            AveWF_09[k, 1:3] = obsweight' * WF_09[ind_k, 1:3]
            y_max = argmax(llhadj_09[k,:])

            BestVals_09[ind_k,:] = [repmat([γ0_δ_vec[y_max,:] llhβ_09[k,y_max]./length(ind_k)],length(ind_k),1) ...
                exp(lip_09[ind_k, y_max]) basellh_09[ind_k,1]...
                γ0_09[ind_k,y_max] γ1_09[ind_k,y_max] γ2_09[ind_k,y_max] ...
                Dm_09[ind_k,y_max] D0_09[ind_k,y_max] pi_09[ind_k,y_max] CSns_09[ind_k,y_max] CSs_09[ind_k,y_max] ]
            RPpost = llhadj12[k,:]' .* vec(imp12)/exp(ltot12[k]-maxtemp12[k])

            ind_k = data12["first"][k]:data12["cdindex"][k]
            WF12[ind_k,1] = pi12[ind_k,:]  * RPpost
            WF12[ind_k,2] = CSns12[ind_k,:]* RPpost
            WF12[ind_k,3] = CSs12[ind_k,:] * RPpost
            obsweight = data12["obsweight"][ind_k,1] ./ sum(data12["obsweight"][ind_k,1])
            AveWF12[k,1:3] = obsweight' * WF12[ind_k, 1:3]

            y_max = argmax(llhadj12[k,:])
            BestVals12[ind_k,:] = [repmat([γ0_δ_vec[y_max,:] llhβ12[k,y_max]./length(ind_k)],length(ind_k),1) ...
                exp(lip12[ind_k,y_max]) basellh12[ind_k,1]...
                γ012[ind_k,y_max] γ112[ind_k,y_max] γ212[ind_k,y_max] ...
                Dm12[ind_k,y_max] D012[ind_k,y_max] pi12[ind_k,y_max] CSns12[ind_k,y_max] CSs12[ind_k,y_max] ]

        end
        fWF.BestVals_09 = [data09["cdid"] data09["numlist"] data09["p"] data09["pdif"] BestVals09]
        fWF.BestVals12 = [data12["cdid"] data12["numlist"] data12["p"] data12["pdif"] BestVals12]
        fWF.BestValsbp = [bp["cdid"] bp["numlist"] bp["p"] repmat([0 0 1 1],bp.N,1) lipb basellhb γ0bp γ1bp γ2bp Dmbp D0bp WFbp]
        fWF.AveWF_09 = AveWF_09
        fWF.AveWF12 = AveWF12
        fWF.WF_09 = WF_09
        fWF.WF12 = WF12
        fWF.WFbp = WFbp
    end

    fother.lip12 = lip12
    fother.llhβ12 = llhβ12
    fother.lip_09 = lip_09
    fother.llhβ_09 = llhβ_09
    fother.γ0_δ_vec = γ0_δ_vec
    fother.imp_09 = imp_09
    fother.imp12 = imp12
    fother.ltot_09 = ltot_09
    fother.ltot12 = ltot12
    fother.imp_09 = imp_09

    fother.γ1_09 = γ1_09
    fother.γ112 = γ112

    fother.lipb = lipb
    fother.γ1bp = γ1bp

    return f, distpara, fother, fWF
end

# Based on file <fullmodelllhWFAug22newtest2015.m>
function full_model(x0, distpara0, gamma0vec, deltavec, data12, data09, bp; WFcal = false)
    # betasigma5 = [gamma0shape gamma0mean09 gamma0mean12 deltasigma
    #               gammaimeanbp betalocal alpha-1 beta gammaishape  gammaimean09
    #               gammaimean12  eta-1 r olp c lamda1 lamda2 betacond betapop
    #               olm oltheta naturaldisappear]
    # betasigma4 = [gamma0shape gamma0mean09 gamma0mean12 deltasigma gammaimeanbp
    #               alpha-1 beta gammaishape  gammaimean09 gammaimean12  eta-1 r
    #               olp c lamda1 lamda2 betacond betapop betalocal olm oltheta]
    betasigma5 = vcat(distpara0, x0[1], x0[2]/(1+x0[1]), x0[3],
        x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1), x0[5]*10*x0[7]/10/8^(-x0[6]-1),
        x0[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1], 0, 0, x0[12:13], x0[14])

    rounderr = 0.025
    naturaldisappear = betasigma5[22]

    betasigma4 = betasigma5[vcat(1:5, 7:19, 6, 20, 21)]

    lamda1   = 0
    lamda2   = 0
    betacond = 0
    betapop  = 0
    olm      = 1.5180
    oltheta  = 8.1041

    lamda1   = betasigma4[15]
    lamda2   = betasigma4[16]
    betacond = betasigma4[17]
    betapop  = betasigma4[18]
    olm      = betasigma4[20]
    oltheta  = betasigma4[21]



    #TODO: ngrid very seldom necessary, needlessly memory intensive
    temp1, temp2   = ndgrid(gamma0vec, deltavec)
    gamma0deltavec = vcat(vec(temp1), vec(temp2))

    Y    = length(temp1)
    M_09 = data09["M"]
    N_09 = data09["N"]
    M_12 = data12["M"]
    N_12 = data12["N"]

    ## Calculation for 09 data
    ltot09    = zeros(M_09)
    llhbeta09 = zeros(M_09, Y)  # record the log likelihood at a fixed beta for a title
    lip09     = zeros(N_09, Y)  # likelihood of each observation at each beta
    gamma109  = lip09
    gamma209  = lip09
    gamma009  = lip09
    D009      = lip09
    Dm09      = lip09
    basellh09 = pdf(Gamma(olm, oltheta), data09["p"]) * 2 * rounderr

    pi09, CSns09, CSs09 = zeros(N_09, Y), zeros(N_09, Y), zeros(N_09, Y)

    # Iterate for gamma0
    @parallel for i = 1:Y
        lip09[:,i], gamma209[:,i], gamma109[:,i], gamma009[:,i], D009[:,i], Dm09[:,i], pi09[:,i], CSns09[:,i], CSs09[:,i] = obscalnewtest2015(
                vcat(gamma0deltavec[i,1], betasigma4[[6 7 8 9 11 12]], lamda1, lamda2, betacond, betapop, 0, betasigma4[13], gamma0deltavec[i,2], betasigma4[14], 1),
                 data09, basellh09, 0, data09["pdif"], rounderr, WFcal)
    end
    lip09, gamma209, gamma109, gamma009, D009, Dm09, pi09, CSns09, CSs09 = obscalnewtest2015([gamma0deltavec[:,1] betasigma4[[6 7 8 9 11 12]] lamda1 lamda2 betacond betapop 0 betasigma4[13] gamma0deltavec[:,2] betasigma4[14] 1],data09,basellh09,0,data09["pdif"],rounderr,WFcal)

    for k=1:M_09
        llhbeta09[k, :] = sum(lip09(data09["first"][k]:data09["cdindex"][k],:))
    end
    maxtemp09 = max(llhbeta09,[],2)
    llhadj09  = exp(llhbeta09 - repmat(maxtemp09,1,Y))

    ## Calculation for 2012 data
    lip12     = zeros(N_12, Y)
    gamma112  = lip12
    gamma212  = lip12
    gamma012  = lip12
    D012      = lip12
    Dm12      = lip12
    ltot12    = zeros(M_12)
    llhbeta12 = zeros(M_12, Y)  # record the log likelihood at a fixed beta for a title
    basellh12 = pdf(Gamma(olm, oltheta), data12["p"]) * 2 * rounderr
    pi12, CSns12, CSs12 = zeros(N_12, Y), zeros(N_12, Y), zeros(N_12, Y)

    # iterate for betas
    @parallel for i = 1:Y
        lip12[:,i], gamma212[:,i], gamma112[:,i], gamma012[:,i], D012[:,i],
        Dm12[:,i], pi12[:,i], CSns12[:,i], CSs12[:,i] = obscalnewtest2015(
            vcat(gamma0deltavec[i,1], betasigma4[[6 7 8 10 11 12]], lamda1, lamda2,
            betacond, betapop, 0, betasigma4[13], gamma0deltavec[i,2], betasigma4[14],
            naturaldisappear), data12, basellh12, 1, data12["pdif"], rounderr, WFcal)
    end
    for k=1:M_12
        llhbeta12[k,:] = sum(lip12(data12["first"][k]:data12["cdindex"][k],:))
    end

    maxtemp12 = max(llhbeta12,[],2)
    llhadj12  = exp(llhbeta12 - fill(maxtemp12,1,Y))

    ## Calculation for 09 offline data
    basellhb = pdf(Gamma(olm, oltheta), bp["p"]) * 2 * rounderr
    function getbmean(gammaimeanbpbetalocal)
        lipb  = obscalnewtest2015([0 betasigma4([6 7 8]) abs(gammaimeanbpbetalocal[1]) betasigma4([ 11 12 ]) lamda1 lamda2 betacond betapop gammaimeanbpbetalocal[2] betasigma4[13] 1 1 1],bp,basellhb,0,bp["p"],rounderr,0)
        return -sum(lipb)
    end

    function integgamma0(gammainput)
        gamma0shape   = gammainput[1]
        gamma0theta09 = gammainput[2]/gammainput[1]
        gamma0theta12 = gammainput[3]/gammainput[1]
        deltamean     = -0.5*gammainput[4]^2
        deltasigma    = gammainput[4]

        # Calculate the importance of the grid points
        gamma0cdfvec = [-gamcdf(gamma0vec[1],gamma0shape,gamma0theta09) gamcdf(gamma0vec,gamma0shape,gamma0theta09) (2 - gamcdf(gamma0vec[end],gamma0shape,gamma0theta09))]
        importance1 = (gamma0cdfvec[3:end]-gamma0cdfvec[1:end-2])/2
        gamma0cdfvec = [-gamcdf(gamma0vec[1],gamma0shape,gamma0theta12) gamcdf(gamma0vec,gamma0shape,gamma0theta12) (2 - gamcdf(gamma0vec[end],gamma0shape,gamma0theta12))]
        importance2 = (gamma0cdfvec[3:end]-gamma0cdfvec[1:end-2])/2
        deltacdfvec = [-normcdf(log(deltavec[1]),deltamean,deltasigma) normcdf(log(deltavec),deltamean,deltasigma) (2 - normcdf(log(deltavec[end]),deltamean,deltasigma))]
        importance3 = (deltacdfvec[3:end]-deltacdfvec[1:end-2])/2
        importance09 = importance1' * importance3
        importance12 = importance2' * importance3

        ltot09 = maxtemp09 + log.(llhadj09 .* vec(importance09))
        ltot12 = maxtemp12 + log.(llhadj12 .* vec(importance12))
        return -(sum(ltot09) + sum(ltot12))
    end

    [distpara1,f1,~,fmindisplay] = fminunc(@integgamma0, betasigma4[1:4], optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    disp(fmindisplay.funcCount)
    [distpara2,f2] = fminunc(@getbmean,betasigma5[5:6],optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    lipb,gamma2bp,gamma1bp,gamma0bp,D0bp,Dmbp,WFbp[:,1],WFbp[:,2],WFbp[:,3] = obscalnewtest2015([0 betasigma4([6 7 8]) abs(distpara2[1]) betasigma4([ 11 12 ]) lamda1 lamda2 betacond betapop distpara2[2] betasigma4[13] 1 1 1],bp,basellhb,0,bp["p"],rounderr,WFcal)

    f = f1 + f2
    distpara = [distpara1 distpara2]

    if WFcal
        WF09       = zeros(N_09, 3)
        WF12       = zeros(N_12, 3)
        AveWF09    = zeros(M_09, 3)
        AveWF12    = zeros(M_12, 3)
        BestVals09 = zeros(N_09,13)
        BestVals12 = zeros(N_12,13)

        for k = 1:M_09
            RPpost = llhadj09[k,:]' .* importance09[:]/exp(ltot09[k] - maxtemp09[k])
            ind_k  = data09["first"][k]:data09["cdindex"][k]

            WF09[ind_k, 1]  = pi09[ind_k,:]*RPpost
            WF09[ind_k, 2]  = CSns09[ind_k,:]*RPpost
            WF09[ind_k, 3]  = CSs09[ind_k,:]*RPpost
            obsweight       = data09["obsweight"][ind_k, 1] ./ sum(data09["obsweight"][ind_k, 1])
            AveWF09[k, 1:3] = obsweight' * WF09[ind_k, 1:3]
            y_max = argmax(llhadj09[k,:])

            BestVals09[ind_k,:] =  [repmat([gamma0deltavec[y_max,:] llhbeta09[k,y_max]./length(ind_k)],length(ind_k),1) ...
                exp(lip09[ind_k, y_max]) basellh09[ind_k,1]...
                gamma009[ind_k,y_max] gamma109[ind_k,y_max] gamma209[ind_k,y_max] ...
                Dm09[ind_k,y_max] D009[ind_k,y_max] pi09[ind_k,y_max] CSns09[ind_k,y_max] CSs09[ind_k,y_max] ]
            RPpost = llhadj12[k,:]' .* vec(importance12)/exp(ltot12[k]-maxtemp12[k])

            ind_k = data12["first"][k]:data12["cdindex"][k]
            WF12[ind_k,1] = pi12[ind_k,:]  * RPpost
            WF12[ind_k,2] = CSns12[ind_k,:]* RPpost
            WF12[ind_k,3] = CSs12[ind_k,:] * RPpost
            obsweight = data12["obsweight"][ind_k,1] ./ sum(data12["obsweight"][ind_k,1])
            AveWF12[k,1:3] = obsweight' * WF12[ind_k, 1:3]

            y_max = argmax(llhadj12[k,:])
            BestVals12[ind_k,:] = [repmat([gamma0deltavec[y_max,:] llhbeta12[k,y_max]./length(ind_k)],length(ind_k),1) ...
                exp(lip12[ind_k,y_max]) basellh12[ind_k,1]...
                gamma012[ind_k,y_max] gamma112[ind_k,y_max] gamma212[ind_k,y_max] ...
                Dm12[ind_k,y_max] D012[ind_k,y_max] pi12[ind_k,y_max] CSns12[ind_k,y_max] CSs12[ind_k,y_max] ]

        end
        fWF.BestVals09 = [data09["cdid"] data09["numlist"] data09["p"] data09["pdif"] BestVals09]
        fWF.BestVals12 = [data12["cdid"] data12["numlist"] data12["p"] data12["pdif"] BestVals12]
        fWF.BestValsbp = [bp["cdid"] bp["numlist"] bp["p"] repmat([0 0 1 1],bp.N,1) lipb basellhb gamma0bp gamma1bp gamma2bp Dmbp D0bp WFbp]
        fWF.AveWF09 = AveWF09
        fWF.AveWF12 = AveWF12
        fWF.WF09 = WF09
        fWF.WF12 = WF12
        fWF.WFbp = WFbp
    end

    fother.lip12 = lip12
    fother.llhbeta12 = llhbeta12
    fother.lip09 = lip09
    fother.llhbeta09 = llhbeta09
    fother.gamma0deltavec = gamma0deltavec
    fother.importance09 = importance09
    fother.importance12 = importance12
    fother.ltot09 = ltot09
    fother.ltot12 = ltot12
    fother.importance09 = importance09

    fother.gamma109 = gamma109
    fother.gamma112 = gamma112

    fother.lipb = lipb
    fother.gamma1bp = gamma1bp

    return f, distpara, fother, fWF
end

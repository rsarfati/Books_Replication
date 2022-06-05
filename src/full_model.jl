# Based on file <fullmodelllhWFAug22newtest2015.m>
function full_model(x0, distpara0, gamma0vec, deltavec, data12, data09, bp; WFcal = false)
    # betasigma5 = [gamma0shape gamma0mean09 gamma0mean12 deltasigma
    # gammaimeanbp betalocal alpha-1 beta gammaishape  gammaimean09 gammaimean12  eta-1 r
    # olp c lamda1 lamda2 betacond betapop olm oltheta naturaldisappear]
    # betasigma4 = [gamma0shape gamma0mean09 gamma0mean12 deltasigma
    # gammaimeanbp  alpha-1 beta gammaishape  gammaimean09 gammaimean12  eta-1 r
    # olp c lamda1 lamda2 betacond betapop betalocal olm oltheta]
    betasigma5 = vcat(distpara0, [x0[1], x0[2]/(1+x0[1]), x0[3],
        x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1), x0[5]*10*x0[7]/10/8^(-x0[6]-1),
        x0[6:11].*[1 0.1 1 0.1 0.01 0.1 ], 0, 0], x0[12:13], x0[14])

    rounderr = 0.025
    naturaldisappear = betasigma5[22]

    betasigma4 = betasigma5[[1:5, 7:19, 6, 20, 21]]

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

    temp1, temp2 = ndgrid(gamma0vec, deltavec)
    gamma0deltavec = [temp1(:) temp2(:)]
    Y = numel(temp1)

    ## Calculation for 09 data
    ltot09    = zeros(data09["M"], 1)
    llhbeta09 = zeros(data09["M"], Y)  # record the log likelihood at a fixed beta for a title
    lip09     = zeros(data09["N"], Y)  # likelihood of each observation at each beta
    gamma109 = lip09
    gamma209 = lip09
    gamma009 = lip09
    D009 = lip09
    Dm09 = lip09
    basellh09 = gampdf(data09["p"], olm, oltheta) * 2 * rounderr
    # basellh09 = gampdf(data09.p,olm,oltheta)
    # p009 = data09.pdif - betacond.*data09.conditiondif
    # p0 = data.p - data.pmin - betacond .* data.conditiondif

    pi09   = zeros(data09.N,Y)
    CSns09 = zeros(data09.N,Y)
    CSs09  = zeros(data09.N,Y)

    # iterate for gamma0

    @parallel for i = 1:Y
        betasigmatemp = betasigma4
        vectemp = gamma0deltavec[i,:]
        [lip09[:,i], gamma209[:,i],gamma109[:,i],gamma009[:,i],D009[:,i],Dm09[:,i],pi09[:,i],CSns09[:,i],CSs09[:,i]] = obscalnewtest2015([vectemp[1] betasigmatemp([6 7 8 9 11 12]) lamda1 lamda2 betacond betapop 0 betasigmatemp[13] vectemp[2] betasigmatemp[14] 1],data09,basellh09,0,data09.pdif,rounderr,WFcal)
        # [mugamma0 alpha-1 beta gammaishape gammaimean eta-1 r lamda1 lamda2 betacond betapop betalocal olp]
    end
    #{
    betasigmatemp = betasigma4
    vectemp = gamma0deltavec[i,:]
    [lip09, gamma209,gamma109,gamma009,D009,Dm09,pi09,CSns09,CSs09] = obscalnewtest2015([gamma0deltavec(:,1) betasigmatemp([6 7 8 9 11 12]) lamda1 lamda2 betacond betapop 0 betasigmatemp[13] gamma0deltavec(:,2) betasigmatemp[14] 1],data09,basellh09,0,data09.pdif,rounderr,WFcal)
    # [mugamma0 alpha-1 beta gammaishape gammaimean eta-1 r lamda1 lamda2 betacond betapop betalocal olp]
    #}


    for k=1:data09.M
        llhbeta09(k,:) = sum(lip09(data09.first(k):data09.cdindex(k),:))
    end
    maxtemp09 = max(llhbeta09,[],2)
    llhadj09 = exp(llhbeta09 - repmat(maxtemp09,1,Y))

    ## Calculation for 12 data

    lip12 = zeros(data12.N,Y)
    gamma112 = lip12
    gamma212 = lip12
    gamma012 = lip12
    D012 = lip12
    Dm12 = lip12
    ltot12=zeros(data12.M,1)
    llhbeta12 = zeros(data12.M,Y)  # record the log likelihood at a fixed beta for a title
    basellh12 =  gampdf(data12.p, olm, oltheta) * 2 * rounderr
    # basellh12 = gampdf(data12.p,olm,oltheta)
    # p012 = data12.pdif  - betacond.*data12.conditiondif

    pi12   = zeros(data12.N,Y)
    CSns12 = zeros(data12.N,Y)
    CSs12  = zeros(data12.N,Y)

    # iterate for betas
    parfor i = 1:Y
        betasigmatemp = betasigma4
        vectemp = gamma0deltavec[i,:]
        [lip12[:,i],gamma212[:,i],gamma112[:,i],gamma012[:,i],D012[:,i],Dm12[:,i],pi12[:,i],CSns12[:,i],CSs12[:,i]] = obscalnewtest2015([vectemp[1] betasigmatemp([6 7 8 10 11 12]) lamda1 lamda2 betacond betapop 0 betasigmatemp[13] vectemp[2] betasigmatemp[14] naturaldisappear],data12,basellh12,1,data12.pdif,rounderr,WFcal)
    end


    for k=1:data12.M
        llhbeta12(k,:) = sum(lip12(data12.first(k):data12.cdindex(k),:))
    end
    maxtemp12 = max(llhbeta12,[],2)
    llhadj12  = exp(llhbeta12 - repmat(maxtemp12,1,Y))

    ## Calculation for 09 offline data

    basellhb = gampdf(bp.p,olm,oltheta)*2*rounderr
    function [ltotb] = getbmean(gammaimeanbpbetalocal)
        #         basellhb= gampdf(bp.p,olm,oltheta)
        lipb = obscalnewtest2015([0 betasigma4([6 7 8]) abs(gammaimeanbpbetalocal[1]) betasigma4([ 11 12 ]) lamda1 lamda2 betacond betapop gammaimeanbpbetalocal[2] betasigma4[13] 1 1 1],bp,basellhb,0,bp.p,rounderr,0)
        ltotb = -sum(lipb)
    end

    function [ltot] = integgamma0(gammainput)

        gamma0shape = gammainput[1]
        gamma0theta09 = gammainput[2]/gammainput[1]
        gamma0theta12 = gammainput[3]/gammainput[1]
        deltamean = -0.5*gammainput[4]^2
        deltasigma = gammainput[4]


        # Calculate the importance of the grid points
        gamma0cdfvec = [-gamcdf(gamma0vec[1],gamma0shape,gamma0theta09) gamcdf(gamma0vec,gamma0shape,gamma0theta09) (2- gamcdf(gamma0vec(end),gamma0shape,gamma0theta09))]
        importance1 = (gamma0cdfvec(3:end)-gamma0cdfvec(1:end-2))/2
        gamma0cdfvec = [-gamcdf(gamma0vec[1],gamma0shape,gamma0theta12) gamcdf(gamma0vec,gamma0shape,gamma0theta12) (2- gamcdf(gamma0vec(end),gamma0shape,gamma0theta12))]
        importance2 = (gamma0cdfvec(3:end)-gamma0cdfvec(1:end-2))/2
        deltacdfvec = [-normcdf(log(deltavec[1]),deltamean,deltasigma) normcdf(log(deltavec),deltamean,deltasigma) (2- normcdf(log(deltavec(end)),deltamean,deltasigma))]
        importance3 = (deltacdfvec(3:end)-deltacdfvec(1:end-2))/2
        importance09 = importance1'*importance3
        importance12 = importance2'*importance3

        ltot09 = maxtemp09 + log(llhadj09*importance09(:))
        ltot09total = -sum(ltot09)
        ltot12 = maxtemp12 + log(llhadj12*importance12(:))
        ltot12total = -sum(ltot12)
        ltot = ltot09total  + ltot12total

    end

    [distpara1,f1,~,fmindisplay] = fminunc(@integgamma0,betasigma4[1:4],optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    disp(fmindisplay.funcCount)
    # [distpara1,f1] = fmincon(@integgamma0,betasigma4[1:4],[],[],[],[],[0 0 0 0],[10 10 10 10],[],optimset('MaxFunEvals',1e4,'Display','off'))
    [distpara2,f2] = fminunc(@getbmean,betasigma5[5:6],optimset('MaxFunEvals',1e4,'Display','off','LargeScale','off'))
    [lipb,gamma2bp,gamma1bp,gamma0bp,D0bp,Dmbp,WFbp(:,1),WFbp(:,2),WFbp(:,3)] = obscalnewtest2015([0 betasigma4([6 7 8]) abs(distpara2[1]) betasigma4([ 11 12 ]) lamda1 lamda2 betacond betapop distpara2[2] betasigma4[13] 1 1 1],bp,basellhb,0,bp.p,rounderr,WFcal)



    #
    # f = f + ltot12total
    # f = f + ltot09total
    f = f1 + f2
    distpara = [distpara1 distpara2]
    #  disp(distpara)
    #
    # integgamma0(distpara)
    # posterior09 = llhadj09.*repmat(importance09(:)',313,1)./repmat(sum(llhadj09.*repmat(importance09(:)',313,1),2),1,3660)
    # posterior12 = llhadj12.*repmat(importance12(:)',313,1)./repmat(sum(llhadj12.*repmat(importance12(:)',313,1),2),1,3660)
    # meanpara09 = posterior09*gamma0deltavec
    # meanpara12 = posterior12*gamma0deltavec
    # bestgamma109 = data09.p.*0
    # bestgamma112 = data12.p.*0
    # [~,bestrandomvec09]=max(llhbeta09,[],2)
    # [~,bestrandomvec12]=max(llhbeta12,[],2)
    # for i = 1:313
    #     bestgamma109(data09.first(i):data09.cdindex(i),1) = gamma109(data09.first(i):data09.cdindex(i),bestrandomvec09(i))
    #     bestgamma112(data12.first(i):data12.cdindex(i),1) = gamma112(data12.first(i):data12.cdindex(i),bestrandomvec12(i))
    # end
    #
    if WFcal

        WF09 = zeros(data09.N,3)
        WF12 = zeros(data12.N,3)
        AveWF09 = zeros(data09.M,3)
        AveWF12 = zeros(data12.M,3)
        BestVals09 = zeros(data09.N,13)
        BestVals12 = zeros(data12.N,13)

        for k = 1:data09.M
            RPpost = llhadj09(k,:)'.*importance09(:)/exp(ltot09(k)-maxtemp09(k))
            Indexk = data09.first(k):data09.cdindex(k)
            WF09(Indexk,1) = pi09(Indexk,:)*RPpost
            WF09(Indexk,2) = CSns09(Indexk,:)*RPpost
            WF09(Indexk,3) = CSs09(Indexk,:)*RPpost
            obsweight = data09.obsweight(Indexk,1)./sum(data09.obsweight(Indexk,1))
            AveWF09(k,1:3) = obsweight'*WF09(Indexk,1:3)
            [~,BestY] = max(llhadj09(k,:))


            BestVals09(Indexk,:) =  [repmat([gamma0deltavec(BestY,:) llhbeta09(k,BestY)./length(Indexk)],length(Indexk),1) ...
                exp(lip09(Indexk,BestY)) basellh09(Indexk,1)...
                gamma009(Indexk,BestY) gamma109(Indexk,BestY) gamma209(Indexk,BestY) ...
                Dm09(Indexk,BestY) D009(Indexk,BestY) pi09(Indexk,BestY) CSns09(Indexk,BestY) CSs09(Indexk,BestY) ]
            RPpost = llhadj12(k,:)'.*importance12(:)/exp(ltot12(k)-maxtemp12(k))
            Indexk = data12.first(k):data12.cdindex(k)
            WF12(Indexk,1) = pi12(Indexk,:)*RPpost
            WF12(Indexk,2) = CSns12(Indexk,:)*RPpost
            WF12(Indexk,3) = CSs12(Indexk,:)*RPpost
            obsweight = data12.obsweight(Indexk,1)./sum(data12.obsweight(Indexk,1))
            AveWF12(k,1:3) = obsweight'*WF12(Indexk,1:3)
            [~,BestY] = max(llhadj12(k,:))
            BestVals12(Indexk,:) = [repmat([gamma0deltavec(BestY,:) llhbeta12(k,BestY)./length(Indexk)],length(Indexk),1) ...
                exp(lip12(Indexk,BestY)) basellh12(Indexk,1)...
                gamma012(Indexk,BestY) gamma112(Indexk,BestY) gamma212(Indexk,BestY) ...
                Dm12(Indexk,BestY) D012(Indexk,BestY) pi12(Indexk,BestY) CSns12(Indexk,BestY) CSs12(Indexk,BestY) ]

            #     obsweight = data09.obsweight(data09.first(k):data09.cdindex(k),1)./sum(data09.obsweight(data09.first(k):data09.cdindex(k),1))
            #     AveWF09(k,1) = ((obsweight'*pi09(data09.first(k):data09.cdindex(k),:)).*llhadj09(k,:))*importance09(:)/exp(ltot09(k)-maxtemp09(k))
            #     AveWF09(k,2) = ((obsweight'*CSns09(data09.first(k):data09.cdindex(k),:)).*llhadj09(k,:))*importance09(:)/exp(ltot09(k)-maxtemp09(k))
            #     AveWF09(k,3) = ((obsweight'*CSs09(data09.first(k):data09.cdindex(k),:)).*llhadj09(k,:))*importance09(:)/exp(ltot09(k)-maxtemp09(k))
            #     obsweight = data12.obsweight(data12.first(k):data12.cdindex(k),1)./sum(data12.obsweight(data12.first(k):data12.cdindex(k),1))
            #     AveWF12(k,1) = ((obsweight'*pi12(data12.first(k):data12.cdindex(k),:)).*llhadj12(k,:))*importance12(:)/exp(ltot12(k)-maxtemp12(k))
            #     AveWF12(k,2) = ((obsweight'*CSns12(data12.first(k):data12.cdindex(k),:)).*llhadj12(k,:))*importance12(:)/exp(ltot12(k)-maxtemp12(k))
            #     AveWF12(k,3) = ((obsweight'*CSs12(data12.first(k):data12.cdindex(k),:)).*llhadj12(k,:))*importance12(:)/exp(ltot12(k)-maxtemp12(k))

        end

        fWF.BestVals09 = [data09.cdid data09.numlist data09.p data09.pdif BestVals09]
        fWF.BestVals12 = [data12.cdid data12.numlist data12.p data12.pdif BestVals12]
        fWF.BestValsbp = [bp.cdid bp.numlist bp.p repmat([0 0 1 1],bp.N,1) lipb basellhb gamma0bp gamma1bp gamma2bp Dmbp D0bp WFbp]
        fWF.AveWF09 = AveWF09
        fWF.AveWF12 = AveWF12
        fWF.WF09 = WF09
        fWF.WF12 = WF12
        fWF.WFbp = WFbp
        # fother.pi09 = pi09
        # fother.CSns09 = CSns09
        # fother.CSs09 = CSs09
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
    # fother.importance12 = importance12
    # fother.posterior09 = posterior09
    # fother.posterior12 = posterior12
    # fother.meanpara09 = meanpara09
    # fother.meanpara12 = meanpara12
    fother.gamma109 = gamma109
    fother.gamma112 = gamma112
    # fother.bestgamma109 = bestgamma109
    # fother.bestgamma112 = bestgamma112

    fother.lipb = lipb
    fother.gamma1bp = gamma1bp

    return f, distpara, fother, fWF
end

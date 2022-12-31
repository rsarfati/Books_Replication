"""
```
full_model(x0::V, distpara0::V, d_on_12::D, d_on_09::D, bp::D;
           γ0vec = vcat(quantile.(Gamma(0.5, 20), 0.005:0.01:0.895), 28:2:60, 64:4:100),
           δ_vec = vcat(exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)), 3:2:20),
           WFcal = false, rounderr = 0.025, parallel = true,
           VERBOSE = true) where {V <: Vector{Float64}, D<:Dict{Symbol, Vector{<:Number}}}
```

Evaluates model at a given x0 and distpara0, given 2009 + 20012 online data (d_on_09, d_on_12),
and offline 2009 data (bp).

Based on the file <fullmodelllhWFAug22newtest2015.m>.
"""
function full_model(x0::V, distpara0::V, d_on_12::D, d_on_09::D, bp::D;
                    γ0vec = [quantile.(Gamma(0.5,20), 0.005:0.01:0.895); 28:2:60; 64:4:100],
                    δ_vec = [exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)); 3:2:20],
                    rounderr = 0.025,
                    # Options
                    WFcal = false, parallel = true,
                    VERBOSE = true) where {V <: Vector{Float64},
                                           D <: Dict{Symbol,Vector{<:Number}}}

    # βσ4 = [γ0shape γ0mean09 γ0mean12 δ_σ γimeanbp
    #        alpha-1 β γishape  γimean09 γimean12  eta-1 r
    #        olp c λ1 λ2 βcond βpop βlocal olm ol_θ]
    βσ4 = [ distpara0[1:5];
            x0[1];
            x0[2]/(1+x0[1]);
            x0[3];
            x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1);
            x0[5]*10*x0[7]/10/8^(-x0[6]-1);
            x0[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1];
            0;
            0;
            distpara0[6];
            x0[12:13]]

    naturaldisappear = x0[14]

    λ1    = βσ4[15]
    λ2    = βσ4[16]
    βcond = βσ4[17]
    βpop  = βσ4[18]
    olm   = βσ4[20]
    ol_θ  = βσ4[21]

    #TODO: ngrid very seldom necessary, needlessly memory intensive
    temp1, temp2 = ndgrid(γ0vec, δ_vec)
    γ0_δ_vec     = hcat(vec(temp1), vec(temp2))

    # Define length constants
    Y = length(temp1)
    N_09, N_12, N_bp = length.([x[:p]       for x in [d_on_09, d_on_12, bp]])
    M_09, M_12, M_bp = length.([x[:d_first] for x in [d_on_09, d_on_12, bp]])

    # Extract variables used with high frequncy
    cdindex_09 = d_on_09[:cdindex]
    first_09   = d_on_09[:d_first]
    p_09       = d_on_09[:p]
    cdindex_12 = d_on_12[:cdindex]
    first_12   = d_on_12[:d_first]
    p_12       = d_on_12[:p]
    disap      = d_on_12[:disappear]

    # Calculation for 09 data
    ltot_09, llhβ_09 = zeros(M_09), zeros(M_09, Y)
    basellh_09 = pdf.(Gamma(olm, ol_θ), p_09) * 2 * rounderr

    ## Calculation for 2012 data
    ltot_12, llhβ_12 = zeros(M_12), zeros(M_12, Y)
    basellh_12 = pdf.(Gamma(olm, ol_θ), p_12) * 2 * rounderr

    ## Calculation for 09 offline data
    basellhb   = pdf.(Gamma(olm, ol_θ), bp[:p]) * 2 * rounderr

    # OPTIMIZE: INIT 2.255119 s (59.73 k alloc: 1.512 GiB, 5.42% gc time) / call
    # WFcal = true:  1.935566 s (53.25 k alloc: 1.267 GiB, 10.41% gc time)
    # WFcal = false: 0.004776 s (33.38 k alloc: 2.988 MiB)
    # WFcal = false: 0.004850 s (33.48 k alloc: 2.991 MiB)
    println(VERBOSE, "Iterating for γ0... (1/2)")
    out_09 = if parallel
       @distributed (hcat) for i in 1:Y
           obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6:9;11:12]]; λ1; λ2; βcond;
                             βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                             d_on_09, N_09, M_09, basellh_09, rounderr;
                             demandcal = false, WFcal = WFcal)
       end
    else
       hcat([obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6:9;11:12]]; λ1; λ2; βcond;
                               βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                               d_on_09, N_09, M_09, basellh_09, rounderr;
                               demandcal = false, WFcal = WFcal) for i=1:Y]...)
    end
    println(VERBOSE, "Completed Iteration for γ0. (2/2)")

    lip_09  = out_09[0*N_09+1:1*N_09,:] # likelihood of each obs at each β
    γ2_09   = out_09[1*N_09+1:2*N_09,:]
    γ1_09   = out_09[2*N_09+1:3*N_09,:]
    γ0_09   = out_09[3*N_09+1:4*N_09,:]
    D0_09   = out_09[4*N_09+1:5*N_09,:]
    Dm_09   = out_09[5*N_09+1:6*N_09,:]
    pi_09   = out_09[6*N_09+1:7*N_09,:]
    CSns_09 = out_09[7*N_09+1:8*N_09,:]
    CSs_09  = out_09[8*N_09+1:9*N_09,:]

    for k = 1:M_09
        # Record log likelihood at a fixed β for a title
        llhβ_09[k,:] .= sum(lip_09[first_09[k] : cdindex_09[k],:])
    end
    maxtemp_09 = maximum(llhβ_09, dims=2)
    llhadj_09  = exp.(llhβ_09 - repeat(maxtemp_09, 1, Y))

    println(VERBOSE, "Iterating for βs... (1/2)")

    out_12 = if parallel
       @distributed (hcat) for i in 1:Y
           obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6:8;10:12]]; λ1; λ2; βcond;
                             βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14];
                             naturaldisappear], d_on_12, N_12, M_12,
                             basellh_12, rounderr;
                             demandcal = true, WFcal = WFcal, disap = disap)
       end
    else
       hcat([obscalnewtest2015([γ0_δ_vec[i,1]; βσ4[[6:8;10:12]]; λ1; λ2; βcond;
                               βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14];
                               naturaldisappear], d_on_12, N_12, M_12,
                               basellh_12, rounderr; demandcal = true,
                               WFcal = WFcal, disap = disap) for i=1:Y]...)
    end
    println(VERBOSE, "Completed Iteration for βs. (2/2)")

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

    getbmean(γ_l) = -sum(obscalnewtest2015([0.; βσ4[6:8]; abs(γ_l[1]);
                                            βσ4[11:12]; λ1; λ2; βcond;
                                            βpop; γ_l[2]; βσ4[13]; 1.; 1.; 1.],
                                            bp, N_bp, M_bp, basellhb, rounderr;
                                            demandcal = false,
                                            WFcal = false)[1:length(basellhb),:])

    function integγ0(γinput::Vector{Float64}; return_all = false)
        γ0shape = γinput[1]
        γ0_θ_09 = γinput[2] / γinput[1]
        γ0_θ_12 = γinput[3] / γinput[1]
        δ_mean  = -0.5 * γinput[4] ^ 2
        δ_σ     = γinput[4]

        # Calculate the importance of the grid points
        γ0_09_cdf_vec = [    -cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec[1]);
                              cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec);
                         2 .- cdf.(Gamma(γ0shape, γ0_θ_09), γ0vec[end])]
        imp1 = (γ0_09_cdf_vec[3:end] - γ0_09_cdf_vec[1:end-2]) ./ 2

        γ0_12_cdf_vec = [    -cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec[1]);
                              cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec);
                         2 .- cdf.(Gamma(γ0shape, γ0_θ_12), γ0vec[end])]
        imp2 = (γ0_12_cdf_vec[3:end] - γ0_12_cdf_vec[1:end-2]) ./ 2

        δ_cdf_vec = [    -cdf( Normal(δ_mean, δ_σ), log( δ_vec[1]));
                          cdf.(Normal(δ_mean, δ_σ), log.(δ_vec));
                     2 .- cdf.(Normal(δ_mean, δ_σ), log( δ_vec[end]))]
        imp3 = (δ_cdf_vec[3:end] - δ_cdf_vec[1:end-2]) ./ 2

        imp_09 = imp1 * imp3'
        imp_12 = imp2 * imp3'

        ltot_09 = maxtemp_09 + log.(llhadj_09 * vec(imp_09) .+ 1e-5)
        ltot_12 = maxtemp_12 + log.(llhadj_12 * vec(imp_12) .+ 1e-5)

        if return_all
            return imp_09, imp_12, ltot_09, ltot_12
        end
        println("ltot_09: ", sum(isnan(ltot_09)))
        println("ltot_12: ", sum(isnan(ltot_12)))
        return -(sum(ltot_09) + sum(ltot_12))
    end

    println(VERBOSE, "Optimizing Pt. I (1/3)")

    res = optimize(integγ0, βσ4[1:4])
    distpara1, f1 = res.minimizer, res.minimum

    println(VERBOSE, "Optimizing Pt. II (2/3)")

    res = optimize(getbmean, distpara0[5:6])
    distpara2, f2 = res.minimizer, res.minimum

    println(VERBOSE, "Finished optimizing! (3/3)")

    out_bp = obscalnewtest2015([0; βσ4[6:8]; abs(distpara2[1]); βσ4[11:12]; λ1; λ2;
                                βcond; βpop; distpara2[2]; βσ4[13]; 1; 1; 1],
                               bp, N_bp, M_bp, basellhb, rounderr;
                               demandcal = false, WFcal = WFcal)

    lipb  = out_bp[0*N_bp+1:1*N_bp,:]
    γ2_bp = out_bp[1*N_bp+1:2*N_bp,:]
    γ1_bp = out_bp[2*N_bp+1:3*N_bp,:]
    γ0_bp = out_bp[3*N_bp+1:4*N_bp,:]
    D0_bp = out_bp[4*N_bp+1:5*N_bp,:]
    Dm_bp = out_bp[5*N_bp+1:6*N_bp,:]

    WF_bp = zeros(N_bp, 3)
    WF_bp[:,1] .= vec(out_bp[6*N_bp+1:7*N_bp,:])
    WF_bp[:,2] .= vec(out_bp[7*N_bp+1:8*N_bp,:])
    WF_bp[:,3] .= vec(out_bp[8*N_bp+1:9*N_bp,:])

    println(VERBOSE, "f1: $f1, f2: $f2")
    f = f1 + f2
    distpara = [distpara1; distpara2]

    fWF    = Dict{String,Any}()
    fother = Dict{String,Any}()

    if WFcal
        println(VERBOSE, "Beginning welfare calculations... (1/2)")

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
            obs_w           = d_on_09[:obs_w][ind_k,1] ./ sum(d_on_09[:obs_w][ind_k,1])
            AveWF_09[k,1:3] = obs_w' * WF_09[ind_k,1:3]
            y_max           = argmax(llhadj_09[k,:])

            BestVals_09[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max, :],
                llhβ_09[k, y_max] ./ length(ind_k)), 1, length(ind_k))',
                exp.(lip_09[ind_k, y_max]), basellh_09[ind_k,1], γ0_09[ind_k,y_max],
                γ1_09[ind_k,y_max], γ2_09[ind_k,y_max],
                Dm_09[ind_k,y_max], D0_09[ind_k,y_max], pi_09[ind_k,y_max],
                CSns_09[ind_k,y_max], CSs_09[ind_k,y_max])

            RPpost = llhadj_12[k,:] .* vec(imp_12) ./ exp(ltot_12[k] - maxtemp_12[k])

            ind_k           = first_12[k]:cdindex_12[k]
            WF_12[ind_k,1]  = pi_12[ind_k,:]  * RPpost
            WF_12[ind_k,2]  = CSns_12[ind_k,:]* RPpost
            WF_12[ind_k,3]  = CSs_12[ind_k,:] * RPpost
            obs_w           = d_on_12[:obs_w][ind_k,1] ./ sum(d_on_12[:obs_w][ind_k,1])
            AveWF_12[k,1:3] = obs_w' * WF_12[ind_k, 1:3]

            y_max = argmax(llhadj_12[k,:])
            BestVals_12[ind_k,:] = hcat(repeat(vcat(γ0_δ_vec[y_max,:],
                llhβ_12[k,y_max] ./ length(ind_k))', length(ind_k), 1),
                exp.(lip_12[ind_k, y_max]), basellh_12[ind_k,1],
                γ0_12[ind_k,y_max], γ1_12[ind_k,y_max], γ2_12[ind_k,y_max],
                Dm_12[ind_k,y_max], D0_12[ind_k,y_max], pi_12[ind_k,y_max],
                CSns_12[ind_k,y_max], CSs_12[ind_k,y_max])
        end

        fWF["BestVals_09"] = hcat(d_on_09[:cdid], d_on_09[:numlist], d_on_09[:p], d_on_09[:pdif], BestVals_09)
        fWF["BestVals_12"] = hcat(d_on_12[:cdid], d_on_12[:numlist], d_on_12[:p], d_on_12[:pdif], BestVals_12)
        fWF["BestValsbp"]  = hcat(bp[:cdid], bp[:numlist], bp[:p],
                                  repeat([0.0, 0.0, 1.0, 1.0]', N_bp, 1),
                                  lipb, basellhb, γ0_bp, γ1_bp, γ2_bp, Dm_bp, D0_bp, WF_bp)
        fWF["AveWF_09"] = AveWF_09
        fWF["AveWF12"]  = AveWF_12
        fWF["WF_09"]    = WF_09
        fWF["WF12"]     = WF_12
        fWF["WFbp"]     = WF_bp

        println(VERBOSE, "Finished welfare calculations. (2/2)")
    end

    fother["lip12"]    = lip_12
    fother["llhβ12"]   = llhβ_12
    fother["imp_12"]   = imp_12
    fother["ltot12"]   = ltot_12
    fother["γ112"]     = γ1_12

    fother["lip_09"]   = lip_09
    fother["llhβ_09"]  = llhβ_09
    fother["imp_09"]   = imp_09
    fother["ltot_09"]  = ltot_09
    fother["γ1_09"]    = γ1_09

    fother["γ0_δ_vec"] = γ0_δ_vec
    fother["lipb"]     = lipb
    fother["γ1bp"]     = γ1_bp

    return f, distpara, fother, fWF, f1, f2
end

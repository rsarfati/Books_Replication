"""
```
full_model(x0::V, distpara0::V, d_on_12::D, d_on_09::D, bp::D;
           γ0vec = [quantile.(Gamma(0.5,20), 0.005:0.01:0.895); 28:2:60; 64:4:100],
           δ_vec = [exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)); 3:2:20],
           rounderr = 0.025, spec = :standard, WFcal = false, parallel = true,
           VERBOSE = true) where {V <: Vector{Float64}, D<:Dict{Symbol, Vector{<:Number}}}
```
Evaluates model at a given x0 and distpara0, given 2009 + 2012 online data (d_on_09, d_on_12),
and offline 2009 data (bp).

Based on the file <fullmodelllhWFAug22newtest2015.m>.
"""
function full_model(x0::V, distpara0::V, d_on_12::D, d_on_09::D, bp::D;
                    γ0vec = [quantile.(Gamma(0.5,20), 0.005:0.01:0.895); 28:2:60; 64:4:100],
                    δ_vec = [exp.(quantile.(Normal(-2,2), 0.01:0.02:0.91)); 3:2:20],
                    rounderr = 0.025,
                    # Options
                    spec = :standard, WFcal = false, parallel = true,
                    VERBOSE = true) where {V <: Vector{Float64},
                                           D <: Dict{Symbol,Vector{<:Number}}}

    # Unpack parameters, perform scaling transformations
    βσ4 = [distpara0[1:5];                          # 1-5: γ0shape γ0mean09 γ0mean12 δ_σ γimeanbp
           x0[1];                                   # 6: alpha-1
           x0[2]/(1+x0[1]);                         # 7: β
           x0[3];                                   # 8: γishape
           x0[4] * 10 * x0[7]/10/9.5^(-x0[6]-1);    # 9: γimean09
           x0[5] * 10 * x0[7]/10/8^(-x0[6]-1);      # 10: γimean12
           x0[6:11] .* [1, 0.1, 1, 0.1, 0.01, 0.1]; # 11-16: eta-1 r olp c λ1 λ2
           0;                                       # 17: βcond
           0;                                       # 18: βpop
           distpara0[6];                            # 19: βlocal
           x0[12:13]]                               # 20-21: olm ol_θ

    nat_disap = x0[14]

    λ1    = βσ4[15]
    λ2    = βσ4[16]
    βcond = βσ4[17]
    βpop  = βσ4[18]
    olm   = βσ4[20]
    ol_θ  = βσ4[21]

    # Parameters for :condition and :cond_list models
    α_c   = (spec != :standard)  ? x0[15] : 0.0 # Condition elasticity (shoppers)
    η_c   = (spec != :standard)  ? x0[16] : 0.0 # Condition coefficient (non-shoppers)
    min_p = (spec == :cond_list) ? x0[17] : 0.0 # Indicator for "lowest-priced listing"

    #OPTIMIZE: ngrid very seldom necessary, needlessly memory intensive
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
    min_p_09   = d_on_09[:has_min_p]
    min_p_12   = d_on_12[:has_min_p]

    # Calculation for 09 online data
    ltot_09, llhβ_09 = zeros(M_09), zeros(M_09, Y)
    basellh_09 = pdf.(Gamma(olm, ol_θ), p_09) * 2 * rounderr

    # Calculation for 2012 online data
    ltot_12, llhβ_12 = zeros(M_12), zeros(M_12, Y)
    basellh_12 = pdf.(Gamma(olm, ol_θ), p_12) * 2 * rounderr

    # Calculation for 09 offline data
    basellhb   = pdf.(Gamma(olm, ol_θ), bp[:p]) * 2 * rounderr

    # OPTIMIZE: INIT 2.255119 s (59.73 k alloc: 1.512 GiB, 5.42% gc time) / call
    # WFcal = true:  1.935566 s (53.25 k alloc: 1.267 GiB, 10.41% gc time)
    # WFcal = false: 0.004776 s (33.38 k alloc: 2.988 MiB)
    # WFcal = false: 0.004850 s (33.48 k alloc: 2.991 MiB)
    println(VERBOSE, "Iterating for γ0... (1/2)")
    out_09 = if parallel
       @distributed (hcat) for i in 1:Y
           obs_cal([γ0_δ_vec[i,1]; βσ4[[6:9;11:12]]; λ1; λ2; βcond;
                    βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                    :d_on_09, N_09, M_09, basellh_09, rounderr;
                    α_c = α_c, η_c = η_c, min_p = min_p,
                    demandcal = false, WFcal = WFcal)
       end
    else
       hcat([obs_cal([γ0_δ_vec[i,1]; βσ4[[6:9;11:12]]; λ1; λ2; βcond;
                      βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; 1],
                     :d_on_09, N_09, M_09, basellh_09, rounderr;
                     α_c = α_c, η_c = η_c, min_p = min_p,
                     demandcal = false, WFcal = WFcal) for i=1:Y]...)
    end
    println(VERBOSE, "Completed Iteration for γ0. (2/2)")

    # Likelihood of each obs at each β
    lip_09 = smooth!(out_09[1:N_09,:])

    # Record log likelihood at a fixed β for a title
    for k = 1:M_09
        llhβ_09[k,:] = sum(lip_09[first_09[k] : cdindex_09[k],:], dims=1)
    end
    maxtemp_09 = maximum(llhβ_09, dims=2)
    llhadj_09  = exp.(llhβ_09 - repeat(maxtemp_09, 1, Y))

    println(VERBOSE, "Iterating for βs... (1/2)")

    out_12 = if parallel
       @distributed (hcat) for i in 1:Y
           obs_cal([γ0_δ_vec[i,1]; βσ4[[6:8;10:12]]; λ1; λ2; βcond;
                    βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; nat_disap],
                    :d_on_12, N_12, M_12, basellh_12, rounderr;
                    α_c = α_c, η_c = η_c, min_p = min_p,
                    demandcal = true, WFcal = WFcal, disap = disap)
       end
    else
        hcat([obs_cal([γ0_δ_vec[i,1]; βσ4[[6:8;10:12]]; λ1; λ2; βcond;
                       βpop; 0; βσ4[13]; γ0_δ_vec[i,2]; βσ4[14]; nat_disap],
                       :d_on_12, N_12, M_12, basellh_12, rounderr;
                       α_c = α_c, η_c = η_c, min_p = min_p,
                       demandcal = true, WFcal = WFcal, disap = disap)
             for i=1:Y]...)
    end
    println(VERBOSE, "Completed Iteration for βs. (2/2)")

    lip_12 = smooth!(out_12[1:N_12,:])

    for k=1:M_12
        llhβ_12[k,:] = sum(lip_12[first_12[k]:cdindex_12[k],:], dims = 1)
    end

    maxtemp_12 = maximum(llhβ_12, dims = 2)
    llhadj_12  = exp.(llhβ_12 .- repeat(maxtemp_12,1,Y))

    getbmean(γ_l) = -sum(obs_cal([0.; βσ4[6:8]; abs(γ_l[1]); βσ4[11:12]; λ1; λ2;
                                  βcond; βpop; γ_l[2]; βσ4[13]; 1.; 1.; 1.],
                                 :bp, N_bp, M_bp, basellhb, rounderr;
                                 α_c = α_c, η_c = η_c, min_p = min_p,
                                 demandcal = false, WFcal = false)[1:length(basellhb),:])

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

        ltot_09 = maxtemp_09 + log.(llhadj_09 * vec(imp_09))
        ltot_12 = maxtemp_12 + log.(llhadj_12 * vec(imp_12))

        if return_all
            return imp_09, imp_12, ltot_09, ltot_12
        end
        return -(sum(ltot_09) + sum(ltot_12))
    end

    println(VERBOSE, "Optimizing Pt. I (1/3)")

    res = optimize(integγ0, 1e-5*ones(4), 50.0*ones(4), βσ4[1:4], Fminbox())
    distpara1, f1 = res.minimizer, res.minimum

    println(VERBOSE, "Optimizing Pt. II (2/3)")

    res = optimize(getbmean, distpara0[5:6])
    distpara2, f2 = res.minimizer, res.minimum

    println(VERBOSE, "Finished optimizing! (3/3)")
    println(VERBOSE, "f1: $f1, f2: $f2")
    println(VERBOSE, "distpara1: $distpara1, distpara2: $distpara2")

    f = f1 + f2
    distpara = [distpara1; distpara2]

    fWF    = Dict{String,Any}()
    fother = Dict{String,Any}()

    if WFcal
        println(VERBOSE, "Beginning welfare calculations... (1/2)")

        # 2009 Online
        CSns_09_full = out_09[1*N_09+1:2*N_09,:]
        CSs_09_full  = out_09[2*N_09+1:3*N_09,:]
        r_p_09_full  = out_09[3*N_09+1:4*N_09,:]
        pi_09_full   = out_09[4*N_09+1:5*N_09,:]

        CSns_09 = Vector{Float64}(undef, N_09)
        CSs_09  = Vector{Float64}(undef, N_09)
        r_p_09  = Vector{Float64}(undef, N_09)
        pi_09   = Vector{Float64}(undef, N_09)
        for k = 1:M_09
            ind_k = first_09[k]:cdindex_09[k]
            y_max = argmax(llhadj_09[k,:])
            CSns_09[ind_k] = CSns_09_full[ind_k, y_max]
            CSs_09[ ind_k] = CSs_09_full[ ind_k, y_max]
            r_p_09[ ind_k] = r_p_09_full[ ind_k, y_max]
            pi_09[  ind_k] = pi_09_full[  ind_k, y_max]
        end

        # 2012 Online
        CSns_12_full = out_12[1*N_12+1:2*N_12,:]
        CSs_12_full  = out_12[2*N_12+1:3*N_12,:]
        r_p_12_full  = out_12[3*N_12+1:4*N_12,:]
        pi_12_full   = out_12[4*N_12+1:5*N_12,:]

        CSns_12 = Vector{Float64}(undef, N_12)
        CSs_12  = Vector{Float64}(undef, N_12)
        r_p_12  = Vector{Float64}(undef, N_12)
        pi_12   = Vector{Float64}(undef, N_12)
        for k = 1:M_12
            ind_k = first_12[k]:cdindex_12[k]
            y_max = argmax(llhadj_12[k,:])
            CSns_12[ind_k] = CSns_12_full[ind_k, y_max]
            CSs_12[ ind_k] = CSs_12_full[ ind_k, y_max]
            r_p_12[ ind_k] = r_p_12_full[ ind_k, y_max]
            pi_12[  ind_k] = pi_12_full[  ind_k, y_max]
        end

        # 2009 Offline
        out_bp = obs_cal([0; βσ4[6:8]; abs(distpara2[1]); βσ4[11:12]; λ1; λ2;
                         βcond; βpop; distpara2[2]; βσ4[13]; 1; 1; 1],
                         :bp, N_bp, M_bp, basellhb, rounderr;
                         α_c = α_c, η_c = η_c, min_p = min_p,
                         demandcal = false, WFcal = WFcal)
         CSns_09_off = out_bp[1*N_bp+1:2*N_bp,:]
         CSs_09_off  = out_bp[2*N_bp+1:3*N_bp,:]
         r_p_09_off  = out_bp[3*N_bp+1:4*N_bp,:]
         pi_09_off   = out_bp[4*N_bp+1:5*N_bp,:]

        if write_output
            @save "$OUTPUT/welfare_estimates_$(string(spec))_$(vint).jld2" CSns_09 CSs_09 r_p_09 pi_09 p_09 CSns_12 CSs_12 r_p_12 pi_12 p_12 CSns_09_off CSs_09_off r_p_09_off pi_09_off
            # 2009 Online
            df_09_on = DataFrame(:title => d_on_09[:cdid], :CSns => vec(CSns_09),
                              :CSs => vec(CSs_09), :r_p => vec(r_p_09),
                              :pi => vec(pi_09), :price => p_09,
                              :condition => d_on_09[:condition],
                              :localint => d_on_09[:localint])
            # 2021 Online
            df_12_on = DataFrame(:title => d_on_12[:cdid], :CSns => vec(CSns_12),
                              :CSs => vec(CSs_12), :r_p => vec(r_p_12),
                              :pi => vec(pi_12), :price => p_12,
                              :condition => d_on_12[:condition],
                              :localint => d_on_12[:localint],
                              :sold => d_on_12[:disappear])
            # 2009 Offline
            df_09_off = DataFrame(:title => bp[:cdid], :CSns => vec(CSns_09_off),
                                :CSs => vec(CSs_09_off), :r_p => vec(r_p_09_off),
                                :pi => vec(pi_09_off), :price => bp[:p],
                                :condition => bp[:condition],
                                :localint => bp[:localint])
            # Write to files
            CSV.write("$OUTPUT/welfare_estimates_09_online_$(string(spec))_$(vint).csv",  df_09_on)
            CSV.write("$OUTPUT/welfare_estimates_12_online_$(string(spec))_$(vint).csv",  df_12_on)
            CSV.write("$OUTPUT/welfare_estimates_09_offline_$(string(spec))_$(vint).csv", df_09_off)
        end
        println(VERBOSE, "Finished welfare calculations. (2/2)")
    end
    return f, x0, distpara, fother, fWF, f1, f2
end

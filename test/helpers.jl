using FileIO, Test, FileIO, JLD2, MAT

path = dirname(@__FILE__) * "/data/"
rng = MersenneTwister(1234)

vecF64(x::Any) = Vector{Float64}(vec(x))
vecI64(x::Any) = Vector{Int64}(vec(x))
vecC64(x::Any) = Vector{ComplexF64}(vec(x))

# Run if you've recently regenerated the .mat files; will reconstruct the .jld2's
reload_matlab_files = true

######################
# demandshopper
######################
if reload_matlab_files
    local v             = matread("$path/testfile_demand.mat")
    local p, cdid, obsw = vecF64.([v["p"], v["data"]["cdid"], v["data"]["obsweight"]])
    local α, β          = Float64.([v["alpha"][1], v["beta"]])

    local f1_mat, f2_mat, f3_mat, sum2_mat, expU_mat = vecF64.([v["f1"],
        v["f2"], v["f3"], v["sum2"], v["expU"]])

    @save "$path/demand_inputs.jld2" α β p cdid obsw
    @save "$path/demand_outputs.jld2" f1_mat f2_mat f3_mat sum2_mat expU_mat
end

@testset "Demand of Shoppers (demandshopper)" begin
    @load "$path/demand_inputs.jld2" α β p cdid obsw
    @load "$path/demand_outputs.jld2" f1_mat f2_mat f3_mat sum2_mat expU_mat

    local f1, f2, f3, sum1, expU = demand_shopper(α, β, p, cdid, obsw; allout = true)

    @test f1 == f1_mat
    @test f2 == f2_mat
    @test f3 == f3_mat
    @test fill(sum1, size(cdid)) == sum2_mat
    @test expU == expU_mat
end

######################
# solve_γ
######################
if reload_matlab_files
    local v    = matread("$path/solvegamma.mat")
    local r    = Float64(v["r"])
    local γ1_m = vecC64(v["gamma1"])
    local Dm, D0, dDm, dD0, γ0, p, l1_m = vecF64.([v["Dm"], v["D0"], v["dDm"], v["dD0"],
                                                   v["gamma0"], v["p"], v["l1"]])
    @save "$path/solvegamma_inputs.jld2" Dm D0 dDm dD0 γ0 p r
    @save "$path/solvegamma_outputs.jld2" γ1_m l1_m
end

@testset "Solve γ" begin
    @load "$path/solvegamma_inputs.jld2" Dm D0 dDm dD0 γ0 p r
    @load "$path/solvegamma_outputs.jld2" γ1_m l1_m

    local γ1, l1 = solve_γ(Dm, D0, dDm, dD0, γ0, p, r; allout = true)

    @test γ1 ≈ γ1_m
    @test l1 == l1_m
end

######################
# welfaresimple
######################
if reload_matlab_files
    local v = matread("$path/testfile_welfare.mat")
    local γ1, γ2, γscale, γ0, olppost, Dm, D0, p0, p, scalarparas,
    pi_mat, CSns_mat, CSs_mat = vecF64.([v["gamma1"], v["gamma2"], v["gammascale"], v["gamma0"],
                                 v["olppost"], v["Dm"], v["D0"], v["pdif"], v["data"]["p"],
                                 v["scalarparas"], v["pi"], v["CSns"], v["CSs"]])
    local cdindex, d_first = vecI64.([v["data"]["cdindex"], v["data"]["first"]])
    local N, M             = Int64.([v["data"]["N"], v["data"]["M"]])

    @save "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
    @save "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat
end

@testset "Welfare Computation (welfaresimple)" begin
    @load "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
    @load "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat

    local pi_test, CSns_test, CSs_test = welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, p0,
                                                       p, N, M, cdindex, d_first, scalarparas)
    @test pi_mat   == pi_test
    @test CSns_mat == CSns_test
    @test maximum(abs.(CSs_mat .- CSs_test)) < 1e-5
end

######################
# welfaresimple
######################
if reload_matlab_files
    v = matread("$path/testfile_welfare.mat")
    local γ1, γ2, γscale, γ0, olppost, Dm, D0, p0, p, scalarparas, pi_mat,
    CSns_mat, CSs_mat = vecF64.([v["gamma1"], v["gamma2"], v["gammascale"], v["gamma0"],
                                 v["olppost"], v["Dm"], v["D0"], v["pdif"], v["data"]["p"],
                                 v["scalarparas"], v["pi"], v["CSns"], v["CSs"]])
    local cdindex, d_first  = Vector{Int64}.(vec.([v["data"]["cdindex"], v["data"]["first"]]))
    local N, M              = Int64.([v["data"]["N"], v["data"]["M"]])

    @save "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
    @save "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat
end

@testset "Welfare Computation (welfaresimple)" begin
    @load "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
    @load "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat

    local pi_test, CSns_test, CSs_test = welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, p0,
                                                 p, N, M, cdindex, d_first, scalarparas)
    @test pi_mat == pi_test
    @test CSns_mat == CSns_test
    @test maximum(abs.(CSs_mat .- CSs_test)) < 1e-5
end

######################
# obscalnewtest2015
######################
if reload_matlab_files
    local v_in = matread("$path/obscal_inputs.mat")
    local βσ3, basellh, p0 = vecF64.([v_in["betasigma3"], v_in["basellh"], v_in["p0"]])
    local data             = v_in["data"]
    local rounderr         = Float64(v_in["rounderr"])

    local v_mat = Dict{String,Any}()
    v_mat["00"] = matread("$path/obscal_outputs_demand=0_WFcal=0.mat")
    v_mat["01"] = matread("$path/obscal_outputs_demand=0_WFcal=1.mat")
    #v_mat["10"] = matread("$path/obscal_outputs_demand=1_WFcal=0.mat")
    #v_mat["11"] = matread("$path/obscal_outputs_demand=1_WFcal=1.mat")

    @save "$path/obscal_inputs.jld2" βσ3 basellh p0 data rounderr
    @save "$path/obscal_outputs.jld2" v_mat
    #@save "$path/obscal_outputs_00.jld2" lip_00 gamma2_00 gamma1_00 gamma0_00 D0_00 Dm_00 pi_00 CSns_00 CSs_00
    #@save "$path/obscal_outputs_01.jld2" lip_01 gamma2_01 gamma1_01 gamma0_01 D0_01 Dm_01 pi_01 CSns_01 CSs_01
    #@save "$path/obscal_outputs_10.jld2" lip_10 gamma2_10 gamma1_10 gamma0_10 D0_10 Dm_10 pi_10 CSns_10 CSs_10
    #@save "$path/obscal_outputs_11.jld2" lip_11 gamma2_11 gamma1_11 gamma0_11 D0_11 Dm_11 pi_11 CSns_11 CSs_11
end

# @load "$path/obscal_outputs_00.jld2" lip_00 gamma2_00 gamma1_00 gamma0_00 D0_00 Dm_00 pi_00 CSns_00 CSs_00
# @load "$path/obscal_outputs_01.jld2" lip_01 gamma2_01 gamma1_01 gamma0_01 D0_01 Dm_01 pi_01 CSns_01 CSs_01
# @load "$path/obscal_outputs_10.jld2" lip_10 gamma2_10 gamma1_10 gamma0_10 D0_10 Dm_10 pi_10 CSns_10 CSs_10
# @load "$path/obscal_outputs_11.jld2" lip_11 gamma2_11 gamma1_11 gamma0_11 D0_11 Dm_11 pi_11 CSns_11 CSs_11

for (i,j) in [(0,0),(0,1)] #,(0,1),(1,1)]
    @testset "obscalnewtest2015 (demandcal = $(i), WFcal = $(j))" begin
    @load "$path/obscal_inputs.jld2" βσ3 basellh p0 data rounderr
    @load "$path/obscal_outputs.jld2" v_mat

    local demandcal = (i == 1)
    local WFcal     = (j == 1)
    local lip, γ2, γ1, γ0, D0, Dm, pi_v, CSns, CSs = obscalnewtest2015(βσ3, data, basellh, p0, rounderr;
                                                         demandcal = demandcal, WFcal = WFcal)

        @show vecF64(v_mat["$i$j"]["lip"])    == lip
        @show vecF64(v_mat["$i$j"]["gamma2"]) == γ2
        @show vecF64(v_mat["$i$j"]["gamma1"]) == γ1
        @show vecF64(v_mat["$i$j"]["gamma0"]) == γ0
        @show vecF64(v_mat["$i$j"]["D0"])     == D0
        @show vecF64(v_mat["$i$j"]["Dm"])     == Dm
        @show vecF64(v_mat["$i$j"]["pi"])     == pi_v
        @show vecF64(v_mat["$i$j"]["CSns"])   == CSns
        @show vecF64(v_mat["$i$j"]["CSs"])    == CSs
    end
end

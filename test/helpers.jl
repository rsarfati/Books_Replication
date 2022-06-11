using FileIO, Test, FileIO, JLD2, MAT

path = dirname(@__FILE__) * "/data/"
rng = MersenneTwister(1234)

# Will reconstruct jld2 files from .mat files
# (Run if you've recently regenerated the .mat files)
reload_matlab_files = true

######################
# demandshopper
######################
if reload_matlab_files
    v = matread("$path/testfile_demand.mat")
    p, cdid, obsw = Vector{Float64}.(vec.([v["p"], v["data"]["cdid"], v["data"]["obsweight"]]))
    α, β = Float64.([v["alpha"][1], v["beta"]])
    f1_mat, f2_mat, f3_mat, sum2_mat, expU_mat = Vector{Float64}.(vec.([v["f1"],
        v["f2"], v["f3"], v["sum2"], v["expU"]]))

    @save "$path/demand_inputs.jld2" α β p cdid obsw
    @save "$path/demand_outputs.jld2" f1_mat f2_mat f3_mat sum2_mat expU_mat
end
@load "$path/demand_inputs.jld2" α β p cdid obsw
@load "$path/demand_outputs.jld2" f1_mat f2_mat f3_mat sum2_mat expU_mat

f1, f2, f3, sum1, expU = demand_shopper(α, β, p, cdid, obsw; testing=true)

@testset "Demand of Shoppers (demandshopper)" begin
    @test f1 == f1_mat
    @test f2 == f2_mat
    @test f3 == f3_mat
    @test fill(sum1, size(cdid)) == sum2_mat
    @test expU == expU_mat
end

######################
# welfaresimple
######################
if reload_matlab_files
    v = matread("$path/testfile_welfare.mat")

    γ1, γ2, γscale, γ0, olppost, Dm, D0, p0, p, scalarparas, pi_mat,
    CSns_mat, CSs_mat = Vector{Float64}.(vec.([v["gamma1"], v["gamma2"], v["gammascale"],
                                             v["gamma0"], v["olppost"], v["Dm"], v["D0"],
                                             v["pdif"], v["data"]["p"], v["scalarparas"],
                                             v["pi"], v["CSns"], v["CSs"]]))
    cdindex, d_first  = Vector{Int64}.(vec.([v["data"]["cdindex"], v["data"]["first"]]))
    N, M              = Int64.([v["data"]["N"], v["data"]["M"]])

    @save "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
    @save "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat
end

@load "$path/welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
@load "$path/welfare_outputs.jld2" pi_mat CSs_mat CSns_mat

pi_test, CSns_test, CSs_test = welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, p0,
                                             p, N, M, cdindex, d_first, scalarparas)

@testset "Welfare Computation (welfaresimple)" begin
    @test pi_mat == pi_test
    @test CSns_mat == CSns_test
    @test maximum(abs.(CSs_mat .- CSs_test)) < 1e-5
end

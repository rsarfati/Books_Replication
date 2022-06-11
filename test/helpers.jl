using FileIO, Test, FileIO, JLD2

path = dirname(@__FILE__) * "/data/"
rng = MersenneTwister(1234)
reload_matlab_files = true

######################
# welfaresimple
######################
if reload_matlab_files
  using MAT
  
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

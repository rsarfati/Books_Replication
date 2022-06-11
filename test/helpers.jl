using FileIO, Test, FileIO, JLD2

rng = MersenneTwister(1234)

reload_matlab_files = true

if reload_matlab_files
  using MAT

  v = matread("/Users/henrygenighx/Desktop/Books/test/testfile_welfare.mat")

  γ1 = Vector{Float64}(vec(v["gamma1"]))
  γ2 = Vector{Float64}(vec(v["gamma2"]))
  γscale = Vector{Float64}(vec(v["gammascale"]))
  γ0 = Vector{Float64}(vec(v["gamma0"]))
  olppost = Vector{Float64}(vec(v["olppost"]))
  Dm = Vector{Float64}(vec(v["Dm"]))
  D0 = Vector{Float64}(vec(v["D0"]))
  p0 = Vector{Float64}(vec(v["pdif"]))
  p = Vector{Float64}(vec(v["data"]["p"]))
  N = Int64(v["data"]["N"])
  M = Int64(v["data"]["M"])
  cdindex = Vector{Int64}(vec(v["data"]["cdindex"]))
  d_first = Vector{Int64}(vec(v["data"]["first"]))
  scalarparas = Vector{Float64}(vec(v["scalarparas"]))

  pi_mat   = Vector{Float64}(vec(v["pi"]))
  CSs_mat  = Vector{Float64}(vec(v["CSs"]))
  CSns_mat = Vector{Float64}(vec(v["CSns"]))

  @save "welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
  @save "welfare_outputs.jld2" pi_mat CSs_mat CSns_mat
end

@load "welfare_inputs.jld2" γ1 γ2 γscale γ0 olppost Dm D0 p0 p N M cdindex d_first scalarparas
@load "welfare_outputs.jld2" pi_mat CSs_mat CSns_mat

pi_test, CSs_test, CSns_test = welfaresimple(γ1, γ2, γscale, γ0, olppost, Dm, D0, p0,
                                             p, N, M, cdindex, d_first, scalarparas)

@test pi_mat == pi_test
@test CSs_mat .≈ CSs_test
@test sum(CSns_mat) ≈ sum(CSns_test)

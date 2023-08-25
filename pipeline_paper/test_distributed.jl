using Pkg
Pkg.activate(".")
using Distributed
using PyCall
using DelimitedFiles

const ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
const ic = readdlm(ic_file)
const r☼ = 8.122
const m = 360.0

# launch worker processes
addprocs(4)

println("Number of processes: ", nprocs())
println("Number of workers: ", nworkers())

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

@everywhere begin
    using PyCall

    pushfirst!(PyVector(pyimport("sys")."path"), "")
    bar = pyimport("test_foo")
    stream = pyimport("stream")
end

# each worker gets its id, process id and hostname
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    θ, ω, β = 40.0, 27.0, 0.002
    param = [θ, ω, β, m, ic, r☼]
    println(id, " " , pid, " ", host, " ", bar.foo(id), stream.chi2_full(θ, ω, β, m, ic, r☼))
end

# remove the workers
for i in workers()
    rmprocs(i)
end
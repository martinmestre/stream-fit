using Pkg
Pkg.activate(".")
using Distributed
using PythonCall
using DelimitedFiles

ic_file = "param_fit_orbit_from_IbataPolysGaiaDR2-data_fixedpot.txt"
icjl = readdlm(ic_file)
r☼ = 8.122
m = 360.0

# launch worker processes
addprocs(3)

println("Number of processes: ", nprocs())
println("Number of workers: ", nworkers())

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

@everywhere begin
    using PythonCall
    pyimport("sys")."path".append("")
    bar = pyimport("test_foo")
    stream = pyimport("stream")
end

# each worker gets its id, process id and hostname
for i in workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    θ, ω, β = 40.0, 27.0, 0.002
    param = [θ, ω, β, m, icjl, r☼]
    χ²s = pyconvert(Float64,bar.foo(θ))
    χ² = pyconvert(Float64,stream.chi2_full(Py(θ), Py(ω), Py(β), Py(m), pyconvert(Vector,icjl), Py(r☼)))
    println(id, " " , pid, " ", host, " ", bar.foo(id), " ", χ²s, " ", χ²)
end

# remove the workers
for i in workers()
    rmprocs(i)
end
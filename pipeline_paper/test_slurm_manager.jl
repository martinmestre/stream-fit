#!/usr/bin/env julia

using Pkg
Pkg.activate(".")
using Distributed, SlurmClusterManager
addprocs(SlurmManager())

@show SlurmManager()
@show nprocs()
@show nworkers()

@everywhere begin
    using Pkg
    Pkg.activate(".")
end

@everywhere println("hello from $(myid()):$(gethostname())")
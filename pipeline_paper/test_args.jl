"""Test arguments"""

println(ARGS)
@show typeof(parse.(Int,ARGS))

using Pkg
ENV["JULIA_GLPK_LIBRARY_PATH"] = "/home/vittorioerba/CNP/glpk-5.0/src/.libs"
Pkg.add(PacckageSpec(url="https://github.com/vittorioerba/CoverFeatureMap.jl"))

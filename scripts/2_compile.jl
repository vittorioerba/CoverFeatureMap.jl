ENV["JULIA_GLPK_LIBRARY_PATH"] = "/home/vittorioerba/CNP/glpk-5.0/src/.libs"
using PackageCompiler, CoverFeatureMap
create_sysimage([:CoverFeatureMap], sysimage_path="./sysimage.so", precompile_execution_file="./3_precompile.jl")


ENV["JULIA_GLPK_LIBRARY_PATH"] = "/home/vittorioerba/CNP/glpk-5.0/src/.libs"
using CoverFeatureMap
cnp(5,5,0.7,true,false,1.,2,1)
cnp(5,5,0.7,true,false,1.,3,1)

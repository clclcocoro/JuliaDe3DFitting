include("JuliaDe3DFitting.jl")

X = [1.2 1 1; 2 2.3 2; 3 3 3.1; 4 4.1 4; 5 5 5]
Y = planeFit(X)
@show X
@show Y

X = [1 0 0; 0 1 0 ; 0 0 1; -1 0 0; 0 -1 0; 0 0 -1]
Y = sphereFit(X, 100)
@show X
@show Y

X = [1 0 0; 0 1 0 ; -1 0 0; 0 -1 0; 1 0 1; 0 1 1; -1 0 1; 0 -1 1; 1 0 2; 0 1 2; -1 0 2; 0 -1 2]
Y = cylinderFit(X)
@show X
@show Y

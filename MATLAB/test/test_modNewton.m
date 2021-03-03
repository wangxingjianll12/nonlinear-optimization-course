clear all; close all; clc; 

% Add path to the functions and Newton algorithm.
addpath('../objective_functions/')
addpath('../algorithms/modNewton/')
addpath('../algorithms/uncMIN/')

%[a, info] = modNewton([1, 0; 0, 2], 1e-5, 1e08);

%params = 0;
funobj = Rosenbrock;
x0 = zeros(funobj.n,1);
[x, info] = uncMIN(funobj, x0)
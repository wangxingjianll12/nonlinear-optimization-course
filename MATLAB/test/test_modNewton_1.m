clear all; close all; clc; 

% Add path to the functions and Newton algorithm.
addpath('../objective_functions/')
addpath('../algorithms/modNewton/')
addpath('../algorithms/uncMIN/')

% Save a dashed line for printing.
dashedline = repelem('-',1,84) ;

% Control parameters in Quadratic Function.
props.n       = 100;
props.density = 0.2;
props.rc      = 1e-2;
props.kind    = 1;
props.g_mean  = 1;
props.g_sd    = 1;

% initial function objects
func_Rosenbrock = Rosenbrock;
func_Genhumps = Genhumps;
func_LeastSquaresTukey_bodyfat = LeastSquaresTukey('../datasets/leastsquares/bodyfat.mat');
func_LeastSquaresTukey_adalone = LeastSquaresTukey('../datasets/leastsquares/abalone.mat');
func_LeastSquaresTukey_bodyfatExpand3 = LeastSquaresTukey('../datasets/leastsquares/bodyfatExpand3.mat');
func_Quadratic = Quadratic(props); 
func_LeastSquares_bodyfat = LeastSquares('../datasets/leastsquares/bodyfat.mat');
func_LeastSquares_abalone = LeastSquares('../datasets/leastsquares/abalone.mat');
func_LeastSquares_bodyfatExpand3 = LeastSquares('../datasets/leastsquares/bodyfatExpand3.mat');
func_Logistic_diabetes = Logistic('../datasets/logistic/diabetes.mat');
func_Logistic_leu = Logistic('../datasets/logistic/leu.mat');
func_Logistic_phishing = Logistic('../datasets/logistic/phishing.mat');  

funobj_names = [
"func_Rosenbrock"
"func_Genhumps"
"func_LeastSquaresTukey_bodyfat"
"func_LeastSquaresTukey_adalone"
"func_LeastSquaresTukey_bodyfatExpand3"
"func_Quadratic"
"func_LeastSquares_bodyfat"
"func_LeastSquares_abalone"
"func_LeastSquares_bodyfatExpand3"
"func_Logistic_diabetes"
"func_Logistic_leu"
"func_Logistic_phishing"
];

funobj = Rosenbrock;
x0 = zeros(funobj.n,1);
params.dir_type   = 'ModifiedNewton';
params.maxit      = 10000;
params.printlevel = 0;
params.tol        = 1e-4;
params.probname   = '';

failed = zeros(0);

for k =1:length(funobj_names)
    funobj_name = funobj_names(k);
    funobj = eval(funobj_name);
    params.probname = char(funobj_name);   % CHANGED for each problem below.
    % Initial estimate of a zero of F.
    x0 = randn(funobj.n, 1);
    % ------------------------------------------
    % Begin testing various objective functions.
    % -------------------------------------------
    fprintf('%s\n',dashedline)
    fprintf(' Begin: Testing of Backtracking Armijo Method\n')
    fprintf('%s\n',dashedline)
    % Call Backtracking Armijo solver.
    [x, info] = uncMIN_1(funobj, x0, params);
    fprintf('Tested %s and exited with status = %2g\n',params.probname, info.status);
    if info.status ~= 0
        failed(k) = 1;
    else 
        failed(k) = 0;
    end
end

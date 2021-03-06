clear;     clc;
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\steepest_descent')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\steepest_problem2')

% Open a file for printing.
outfileID = fopen('test_Project2_2.out','w+');

% Control parameters in structure params (Used for all problems)
params.maxiter    = 1e+5;  
params.printlevel = 1;      
params.tol        = 1e-4;
params.alphainit  = 1e-5; 
params.tau        = 0.8;         % tau in (0,1)
params.eta        = 1.2;         % eta > 1
params.outfileID  = outfileID; 

% Determine the objective function: 1   - Rosenbrock
%                                   2   - Genhumps
%                                   3   - Quadratic
%                                   4.1 - Least-Squares-Tukey bodyfat
%                                   4.2 - Least-Squares-Tukey abalone
%                                   4.3 - Least-Squares-Tukey bodyfatExpand3
%                                   5.1 - Least-Squares bodyfat
%                                   5.2 - Least-Squares abalone
%                                   5.3 - Least-Squares bodyfatExpand3
%                                   6.1 - Logistic Regression diabetes
%                                   6.2 - Logistic Regression leu
%                                   6.3 - Logistic Regression phishing
det = 6.3;

% Save a dashed line for printing.
dashedline = repelem('-',1,93);

% ------------------------------------------
% Begin testing various objective functions.
% -------------------------------------------
fprintf('%s\n',dashedline)
fprintf(' Begin: Testing of steepest descent solver\n')
fprintf('%s\n',dashedline)

switch det
case 1
% ------------------------------------------
% Test function: Rosenbrock function
% ------------------------------------------

% Gather the object Rosenbrock.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Rosenbrock............')
funobj = Rosenbrock;

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minimizer of f.
x0        = [4;6];

% Control parameters in structure params.
params.probname   = 'Rosenbrock'; % CHANGES for each new problem below.

% Call steepest descent solver for problem 2.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 2
% ------------------------------------------
% Test: Genhumps function
% ------------------------------------------
    
% Gather the object Genhumps.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Genhumps..............')
funobj = Genhumps;

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a zero of F.
x0 = ones(5,1);

% Name of function.
params.probname = 'Genhumps';

% Call steepest descent solver for problem 2.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 3
% ------------------------------------------
% Test: Quadratic function
% ------------------------------------------

% Gather the object Quadratic.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Quadratic.............')
props.n       = 100;
props.density = 0.2;
props.rc      = 1e-2;
props.kind    = 1;
props.g_mean  = 1;
props.g_sd    = 1;
funobj        = Quadratic(props);

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a zero of F.
x0 = ones(funobj.n,1);

% Name of function.
params.probname = 'Quadratic (100,0.2,1e-3,1,1,1)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 4.1
% ------------------------------------------
% Test function: Least-Squares-Tukey
% ------------------------------------------

% bodyfat data set
% -----------------

% Gather the object Least-Squares-Tukey.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares-Tukey...')
funobj = LeastSquaresTukey('../datasets/leastsquares/bodyfat.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares-Tukey (data:bodyfat)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 4.2
% ------------------------------------------
% Test function: Least-Squares-Tukey
% ------------------------------------------

% abalone data set
% -----------------

% Gather the object Least-Squares-Tukey.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares-Tukey...')
funobj = LeastSquares('../datasets/leastsquares/abalone.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares-Tukey (data:abalone)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 4.3
% ------------------------------------------
% Test function: Least-Squares-Tukey
% ------------------------------------------

% bodyfatExpand3 data set
% ------------------------

% Gather the object Least-Squares-Tukey.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares-Tukey...')
funobj = LeastSquares('../datasets/leastsquares/bodyfatExpand3.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares-Tukey (data:bodyfatExpand3)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 5.1
% ------------------------------------------
% Test function: Least-Squares
% ------------------------------------------

% bodyfat data set
% -----------------

% Gather the object Least-Squares.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares.........')
funobj = LeastSquares('../datasets/leastsquares/bodyfat.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares (data:bodyfat)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 5.2
% ------------------------------------------
% Test function: Least-Squares
% ------------------------------------------

% abalone data set
% -----------------

% Gather the object Least-Squares.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares.........')
funobj = LeastSquares('../datasets/leastsquares/abalone.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares (data:abalone)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 5.3
% ------------------------------------------
% Test function: Least-Squares
% ------------------------------------------

% bodyfatExpand3 data set
% ------------------------

% Gather the object Least-Squares.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Least-Squares.........')
funobj = LeastSquares('../datasets/leastsquares/bodyfatExpand3.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minmizer of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Least-Squares (data:bodyfatExpand3)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 6.1
% ------------------------------------------
% Test function: Logistic Regression
% ------------------------------------------

% diabetes data set
%-------------------

% Gather the object Logistic
fprintf(' Testing algorithm STEEPEST_DESCENT on function Logistic..............')
funobj = Logistic('../datasets/logistic/diabetes.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minimize of f.
x0 = ones(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Logistic (data:diabetes)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 6.2
% ------------------------------------------
% Test function: Logistic Regression
% ------------------------------------------

% leu data set
%-------------------

% Gather the object Logistic
fprintf(' Testing algorithm STEEPEST_DESCENT on function Logistic..............')
funobj = Logistic('../datasets/logistic/leu.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minimize of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Logistic (data:leu)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

case 6.3
% ------------------------------------------
% Test function: Logistic Regression
% ------------------------------------------

% phishing data set
%-------------------

% Gather the object Logistic
fprintf(' Testing algorithm STEEPEST_DESCENT on function Logistic..............')
funobj = Logistic('../datasets/logistic/phishing.mat');

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Initial estimate of a minimize of f.
x0 = zeros(size(funobj.A,2), 1);

% Name of function.
params.probname = 'Logistic (data:phishing)';

% Call steepest descent solver.
[~,info] = steepest_problem2(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info.status);

end

% ------------
% Finish up.
% ------------
fprintf('%s\n',dashedline)
fprintf(' End: Testing of STEEPEST_DESCENT solver\n')
fprintf('%s\n',dashedline)
fclose('all');
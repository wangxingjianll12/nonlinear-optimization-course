% clear;     clc;
det = 1.2;
switch det
case 1.1
% Add path to the functions and Newton algorithm.
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\newton')

% Save a dashed line for printing.
dashedline = repelem('-',1,84) ;

% ------------------------------------------
% Begin testing various objective functions.
% -------------------------------------------
fprintf('%s\n',dashedline)
fprintf(' Begin: Testing of Newton solver\n')
fprintf('%s\n',dashedline)

% ------------------------------------------
% Test: Rosenbrock function
% ------------------------------------------

% Gather the object Rosenbrock.
fprintf(' Testing algorithm NEWTON on function Rosenbrock.............')
funobj = Rosenbrock;

% Define function handles for computing F and its Jacobian J.
Ffunc = @funobj.grad;
Jfunc = @funobj.hess;

% Initial estimate of a zero of F.
x0 = [150;200];

% Open a file for printing.
outfileID = fopen('test_newton.out','w+');

% Control parameters in structure params.
params.maxiter    = 30;            % Used for all problems.
params.printlevel = 1;             % Used for all problems.
params.tol        = 1e-5;          % Used for all problems.
params.outfileID  = outfileID;     % Used for all problems.
params.probname   = 'Rosenbrock';  % CHANGED for each problem below.

% Call Newton Method solver.
[x1,info1] = newton(Ffunc,Jfunc,x0,params);
fprintf('exited with status = %2g\n',info1.status);

case 1.2
% Add path to the functions and Newton algorithm.
addpath('../objective_functions/')
addpath('../algorithms/steepest_descent/')

% Save a dashed line for printing.
dashedline = repelem('-',1,93);

% ------------------------------------------
% Begin testing various objective functions.
% -------------------------------------------
fprintf('%s\n',dashedline)
fprintf(' Begin: Testing of steepest descent solver\n')
fprintf('%s\n',dashedline)

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
x0 = [50;100];

% Open a file for printing.
outfileID = fopen('test_steepest_descent.out','w+');

% Control parameters in structure params.
params.maxiter    = 1e+5;         % Used for all problems.
params.printlevel = 1;            % Used for all problems.
params.tol        = 1e-5;         % Used for all problems.
params.stepchoice = 'fixed';      % Used for all problems.
params.stepsize   = 1e-6;         % Used for all problems.
params.outfileID  = outfileID;    % Used for all problems.
params.probname   = 'Rosenbrock'; % CHANGES for each new problem below.

% Call steepest descent solver.
[x2,info2] = steepest_descent(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info2.status);

case 2
% Add path to the functions and Steepest Descent method.
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\steepest_descent')

step = [5e-5 1e-4 3e-4 5e-4 8e-4 1e-3 1.2e-3 1.5e-3 1.7e-3 1.8e-3 1.9e-3 2e-3];
x_rec    = [0;0];
iter_rec = 0;
f_rec    = 0;

for i = 1:length(step)
% Save a dashed line for printing.
dashedline = repelem('-',1,93);

% ------------------------------------------
% Begin testing various objective functions.
% -------------------------------------------
fprintf('%s\n',dashedline)
fprintf(' Begin: Testing of steepest descent solver\n')
fprintf('%s\n',dashedline)

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
x0 = [1;2];

% Open a file for printing.
outfileID = fopen('test_steepest_descent.out','w+');

% Control parameters in structure params.
params.maxiter    = 1e+6;         % Used for all problems.
params.printlevel = 1;            % Used for all problems.
params.tol        = 1e-5;         % Used for all problems.
params.stepchoice = 'fixed';      % Used for all problems.
params.stepsize   = step(i);      % Used for all problems.
params.outfileID  = outfileID;    % Used for all problems.
params.probname   = 'Rosenbrock'; % CHANGES for each new problem below.

% Call steepest descent solver.
[x2,info2] = steepest_descent(f_hand,g_hand,x0,params);
fprintf('exited with status = %2g\n',info2.status);

% Record the useful information
x_rec(:,i)  = x2;
iter_rec(i) = info2.iter;
f_rec(i)    = info2.f;

end

figure(1)
subplot(2,1,1);
semilogx(step,f_rec,'*-');
ylim([2.46e-5 2.55e-5]);
xlabel('Step size \alpha');
ylabel('norm(f)');
subplot(2,1,2);
loglog(step,iter_rec,'*-');
ylim([1e3 1e6]);
xlabel('Step size \alpha');
ylabel('Number of iteration');

case 3
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\steepest_descent')

% ------------------------------------------
% Test function: Rosenbrock function
% ------------------------------------------

% Gather the object Rosenbrock.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Rosenbrock............')
funobj = Rosenbrock;

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

% Range of coordinates
x1 = -1:0.1:1;
x2 = -10:1:10;
info1 = 0;        % record the status
info2 = 0;        % record the max alpha

for i = 1:length(x1)
for j = 1:length(x2)
% Initial estimate of a minimizer of f.
x0    = [x1(i);x2(j)];
iter  = 0;

x        = x0;
xprev    = x;
norm_x   = norm(x);
f        = feval(f_hand,x);
g        = feval(g_hand,x);
tol      = 2;
maxiter  = 1e6;
stepsize = 1e-8;

while (1)
    if iter > maxiter
        status = 2;          % max iteration reached
        info1(i,j) = status;
        info2(i,j) = alpha;
        break
    end
    if iter > 1
        diffrec = fprev-alpha/2*norm(gprev)^2-f;
    end
    x     = xprev;
    fprev = feval(f_hand,x);
    gprev = feval(g_hand,x);
    
    iter  = iter+1;
    alpha = iter*stepsize;
    p     = -alpha*gprev;
    x     = x+p;
    f     = feval(f_hand,x);
    g     = feval(g_hand,x);
    diff = fprev-alpha/2*norm(gprev)^2-f;
    
    % record the first difference as the convergence
    if iter == 1
        diff0 = diff;
    end
    
    if f > fprev-alpha/2*norm(gprev)^2
        if -diff > tol*diff0
            status = 1;       % stepsize too large
            info1(i,j) = status;
            info2(i,j) = alpha-stepsize;
            break
        else
            status = 0;       % good result
            info1(i,j) = status;
            info2(i,j) = alpha-stepsize;
            break
        end
    end
end

fprintf('exited with status = %2g\n',status);

end
end

case 4
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\objective_functions')
addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\steepest_descent')

% ------------------------------------------
% Test function: Rosenbrock function
% ------------------------------------------

% Gather the object Rosenbrock.
fprintf(' Testing algorithm STEEPEST_DESCENT on function Rosenbrock............')
funobj = Rosenbrock;

% Define function handles for computing F and its Jacobian J.
f_hand = @funobj.func;
g_hand = @funobj.grad;

x0        = [4;6];
iter      = 0;
alphainit = 1e-2;

x        = x0;
norm_x   = norm(x);
maxiter  = 1e6;
norm_x   = norm(x);
f        = feval(f_hand,x);
g        = feval(g_hand,x);
norm_g   = norm(g);
norm_g0  = norm_g; % Save value at x0 to use in relative stopping condition.
tol      = 1e-4;
tau      = 0.8;
eta      = 1.2;

while (1)
    iter = iter+1;
% Check for termination
    if ( norm_g <= tol*max(1,norm_g0) )
       status = 0;
       outcome = ' Relative stopping tolerance reached';
       break
    elseif ( iter >= maxiter )
       status = 2;
       outcome = ' Maximum allowed iterations reached';
       break
    end

    if iter == 1
        alpha = alphainit;
        p     = -alpha*g;
        x_new = x+p;
        diff  = feval(f_hand,x_new)-f+alpha/2*norm(g)^2;
        if diff > 0
            while (1)
                alpha = alpha*tau;
                p     = -alpha*g;
                x_new = x+p;
                diff  = feval(f_hand,x_new)-f+alpha/2*norm(g)^2;
                if diff <= 0
                    break
                end
            end
        else
            while (1)
                p     = -alpha*eta*g;
                x_new = x+p;
                diff  = feval(f_hand,x_new)-f+alpha/2*norm(g)^2;
                if diff > 0
                    break
                end
                alpha = alpha*eta;
            end
        end
        p      = -alpha*g;
        norm_p = alpha*norm_g;
        x      = x + p;
        norm_x = norm(x);
        f      = feval(f_hand,x);
        g      = feval(g_hand,x);
        norm_g = norm(g);
    else
        p     = -alpha*g;
        x_new = x + p;
        diff  = feval(f_hand,x_new)-f+alpha/2*norm(g)^2;
        if diff > 0
            while (1)
                alpha = alpha*tau;
                p     = -alpha*g;
                x_new = x+p;
                diff  = feval(f_hand,x_new)-f+alpha/2*norm(g)^2;
                if diff <= 0
                    break
                end
            end
        end
        p      = -alpha*g;
        norm_p = alpha*norm_g;
        x      = x + p;
        norm_x = norm(x);
        f      = feval(f_hand,x);
        g      = feval(g_hand,x);
        norm_g = norm(g);
    end
    
% Check for NaNs in F and J.
    if isnan(f) || isinf(f)
        status = -2;
        outcome = ' ERROR (NaN/Inf when evaluating f)';
        break
    elseif sum(isnan(g(:))) >= 1 || sum(isinf(g(:))) >= 1
        status = -2;
        outcome = ' ERROR (NaN/Inf when evaluating g)';
        break
    end
end

fprintf('exited with status = %2g\n',status);

end
function [x, info] = uncMIN_1(fun_hands, x0, params)
%------------------------------------------------------------------------
% The function call
%      [x,info] = uncMIN(fun_hands, x0, params )
% aims to compute a vector x that solves the optimization problem
%      minimize f(x)
% usign Backtracking Armijo Method.  The arguments are defined as follows:
%
%  Input arguments:
%  ----------------
%    f_hand     : structure with the following function members:
%                 func       : handle of a function that computes f(x).
%                 grad       : handle of a function that computes gradient
%                              of f(x) 
%                 hess       : handle of a function that computes hessian
%                              of f(x)
%                 hessvecprod: handle of a function that computes the
%                              hessian of f(x) and any vector product
%    x0         : initial guess of a minimizer of f.
%    params     : structure with the following members:
%                 dir_type   : method used to calculate the descent
%                              direction. It can either take the value
%                              SteepestDescent or ModifiedNewton
%                 maxit      : maximum number of allowed iterations.
%                 printlevel : amount of printing to perform.
%                                 0  no printing
%                                 1  single line of output per iteration
%                              When printing, the following is displayed:
%                                 Iter       current iteration number
%                                 Norm-F     2-norm of F at current x
%                                 Norm-x     2-norm of current iterate x
%                                 Norm-step  2-norm of Newton step
%                 tol        : desired stoppping tolerance.
%                 probname   : problem name used for printing purposes.
%
%  Output arguments:
%  -----------------
%    x      : final iterate computed
%    info   : structure with the following members:
%             f         : value of f at the final iterate.
%             f_evals   : the total number of function evaluations.
%             g         : value of the gradient of f at the final iterate.
%             g_evals   : the total number of gradient evaluations.
%             g_norm    : norm of the gradient of f at the final iterate.
%             H         : the hessian of function evaluated at the final
%                         iterate if the algorithm uses one
%             H_evals   : the total number of hessian evaluations
%             HV_evals  : the total number of hessian vector product evaluations
%             iter      : total number of iterations performed.
%             status    : integer indicating the reason for termination
%                         0  Successful termination. A problem is considered 
%                            to have been successful solved if it finds an x
%                            satisfying ||g|| <= tol*max(1,||g0)||) whgere g
%                            is as above and g0 is the gradient of f evaluated
%                            at the input parameter x0.
%                         1  Step computed is too small to make more progress
%                            or maximum number of iterations is reached or an
%                            error in the inputs is detected or a NAN or Inf
%                            is detected or a warning is encountered

addpath('../modNewton')

% Set dummy values for outputs; prevents errors resulting from bad inputs.  
x       = [];
f       = [];  info.f = f;
f_evals = 0 ;  info.f_evals = f_evals;
g       = [];  info.g = g;
g_evals = 0 ;  info.g_evals = g_evals;  
g_norm  = 0 ;  info.g_norm = g_norm;
H       = [];  info.H = H;
H_evals = 0 ;  info.H_evals = H_evals;
Hv_evals= 0 ;  info.Hv_evals = Hv_evals; 
iter    = 0 ;  info.iter = iter;
iter_alp= zeros(0);  info.iter_alp = iter_alp;
status  = 1 ;  info.status = status;

%Check whether enough input parameters were supplied 
if nargin < 2
    info.status = 1;
    fprintf('\n Backtracking Armijo Linesearch: not enough input arguments \n')
    return
end

% If two input arguments, then define default values for third input.
if nargin == 2
    params.dir_type   = 'ModifiedNewton';
    params.maxit      = 10000;
    params.printlevel = 1;
    params.tol        = 1e-8;
    params.probname   = '';
end

% Check that the initial point has the desired length.
if length(x0) ~= fun_hands.n
    fprintf('\n Check the length of input argument x0 \n');
    info.status = 1;
    return
end

% Check all control parameters
if isfield(params, 'dir_type')
    method = params.dir_type;
    if strcmp(method, 'SteepestDescent') ~= 1 && strcmp(method, 'ModifiedNewton') ~=1
        fprintf('\n Direction selection method should be valid \n');
        info.status = 1;
        return
    end
else 
    fprintf('\n Control parameter dir_type is not supplied \n')
    info.status = 1;
    return
end

if isfield(params, 'maxit')
    parameter = params.maxit;
    if parameter <= 0
        fprintf('\n Maximum iterations have to greater than zero \n')
        info.status = 1;
        return
    end
else
    fprintf('\n Control parameter maxit is not supplied \n')
    info.status = 1;
    return
end

if isfield(params, 'printlevel')
    parameter = params.printlevel;
    if parameter < 0 || parameter > 1
        fprintf('\n Print level can either be 0 or 1 \n')
        info.status = 1;
        return
    end
else
    fprintf('\n Control parameter printlevel is not supplied \n')
    info.status = 1;
    return
end

if isfield(params, 'tol')
    parameter = params.tol;
    if parameter < 0 
        fprintf('\n Tolerance has to be greater than zero \n')
        info.status = 1;
        return
    end
else
    fprintf('\n Control parameter tol is not supplied \n')
    info.status = 1;
    return
end
    
if isfield(params,'probname')
    parameter = params.probname;
    if ~ischar(parameter)
        fprintf('\n Invalid control parameter probname. \n');
        info.status = 1;
        return
    end 
end

% Initialize the parameters
x       = x0;
f       = fun_hands.func(x0);    
g       = fun_hands.grad(x0);      
g_norm  = norm(g);  
H       = fun_hands.hess(x0);  

% Backtracking Armijo parameters
alpha = 1; % Initial stepsize
tau = 0.5; % Backtracking parameter
nu = 0.0001; % Bounded (0, 1)

% Modified Newtons parameters
lam_min = 1e-05;
lam_max = 1e08;

% Save things to make printing easier.
dashedline = repelem('-',1,81) ;
header     = '  Iter         f          norm_g         norm_x  ';

% Print column header and value of f at initial point.
if params.printlevel ~= 0   
    fprintf('%s\n',dashedline);
    fprintf('                            Backtracking Armijo Method \n');
    fprintf('%s\n',dashedline);
    fprintf(' maximum iterations            : %g\n',params.maxit);
    fprintf(' print level                   : %g\n',params.printlevel);
    fprintf(' termination tolerance         : %1.2e\n',params.tol);
    fprintf(' problem name                  : %s\n',params.probname);
    fprintf(' Direction methos used         : %s\n',params.dir_type);
    fprintf('%s\n',dashedline);
    fprintf('%s\n',header);
    fprintf(' %5g %14.7e %14.7e %14.7e', iter, f, g_norm, norm(x));
end

% Constant
TINY = eps^(2/3); % Determines if step is too small to make progress.

% Run the Backtrack Armijo Algorithm with either steepest descent or
% modified newtons
for i = 1:params.maxit+10
    
    % Check stopping conditions
    if g_norm <= params.tol*max(1, fun_hands.grad(x0))
        fprintf('\n Maximum tolerance reached. \n')
        info.status = 1;
        break
    elseif iter == params.maxit
        fprintf('\n Maximum iterations reached. \n')
        info.status = 1;
        break
    end
    
    % Choose descent direction (only two available)
    if strcmp(params.dir_type,'SteepestDescent') 
        p = - g;
    else
        p = -modNewton(H, lam_min, lam_max)\g;
        H_evals = H_evals + 1;
    end
    
    %Set initial guess for alpha and iterations
    alpha_new = alpha;
    
    % Choose stepsize
    for j = 1:10*params.maxit
        % Choose next point
        x_new = x + alpha_new*p;
        f_new = fun_hands.func(x_new);
        g_new = fun_hands.grad(x_new);
        
        if f_new <= f + nu*alpha_new*g_new'*p
            alpha = alpha_new;
            break;
        else  
            alpha_new = alpha_new*tau;
        end
        iter_alp(i) = j;
    end
    
    norm_p = alpha*norm(p);
    norm_xprev = norm(x);

    %Take next step
    x = x + alpha*p;
    f = fun_hands.func(x);
    g = fun_hands.grad(x);      
    g_norm = norm(g);  
    iter = iter + 1;
    f_evals = f_evals +1;
    g_evals = g_evals + 1;

    % Print iterate information, if needed.
    if params.printlevel ~= 0
        fprintf(' %14.7e \n', norm_p); % End previous line.
        fprintf(' %5g %14.7e %14.7e %14.7e', iter, f, g_norm, norm(x)); % Start next line.
    end
    
    % Check for NaNs in F and g. Check for zero values of gradient
    if isnan(f) || isinf(f)
        info.status = 1;
        fprintf('\n ERROR (NaN/Inf when evaluating f) \n')
        break
    elseif sum(isnan(g(:))) >= 1 || sum(isinf(g(:))) >= 1
        info.status = 1;
        fprintf('\n ERROR (NaN/Inf when evaluating g) \n');
        break
    elseif g_norm == 0
        info.status = 0;
        fprintf('\n The norm of gradient is zero. \n')
        return
    end

    if norm_p < (1 + norm_xprev)*TINY
        info.status = 1;
        fprintf('\n  Step is too small to make additional progress. \n')
        break
    end
end

% Fill output variable inform.
info.f = f;
info.f_evals = f_evals;
info.g = g;
info.g_evals = g_evals;  
info.g_norm = g_norm;
info.H = H;
info.H_evals = H_evals;
info.Hv_evals = Hv_evals; 
info.iter = iter;
info.iter_alp = iter_alp;

end
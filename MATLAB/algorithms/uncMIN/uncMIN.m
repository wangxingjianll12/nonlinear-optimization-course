function [x,info] = uncMIN(fun_hands,x0,params)
%------------------------------------------------------------------------
% The function call
%      [x,info] = uncMIN(fun_hands,x0,params)
% aims to compute a vector x that solves the optimization problem
%      minimize f(x)
% using the backtracking-Armijo linesearch method.  The arguments are 
% defined as follows:
%
%  Input arguments:
%  ----------------
%    fun_hands  : structure with members f_hand, g_hand, H_hand and Hv_hand
%                 that are all function handles.
%    x0         : initial guess of a minimizer of f.
%    params     : structure with the following members:
%                 dir_type   : a string that indicates the type of search
%                              direction that should be computed. Two 
%                              values are allowed, the SteepestDescent and 
%                              ModifiedNewton.
%                 maxiter    : maximum number of allowed iterations.
%                 printlevel : amount of printing to perform.
%                                 0  no printing
%                                 1  single line of output per iteration
%                              When printing, the following is displayed:
%                                 Iter       current iteration number
%                                 Norm-F     2-norm of F at current x
%                                 Norm-x     2-norm of current iterate x
%                                 Norm-step  2-norm of Newton step
%                 tol        : desired stoppping tolerance.
%                 outputfile : name of a file for output to be printed.
%                 probname   : problem name used for printing purposes.
%  Output arguments:
%  -----------------
%    x     : final iterate computed
%    info  : structure with the following members:
%            f        : value of f at the final iterate
%            f_evals  : the total number of function evaluations
%            g        : value of the gradient of f at the final iterate
%            g_evals  : the total number of gradient evaluations
%            H        : value of the Hessian of f at the final iterate
%            H_evals  : the totalnumber of Hessian evaluations
%            Hv_evals : the totalnumber of Hessian-vector products
%            iter     : total number of iterations performed
%            status   : integer indicating the reason for termination
%                       0  Successful termination. A problem is considered 
%                          to have been successful solved if it finds an x
%                          satisfying ||g|| <= tol*max(1,||g0)||) whgere g
%                          is as above and g0 is the gradient of f evaluated
%                          at the input parameter x0.
%                       1  Step computed is too small to make more progress
%                          or maximum number of iterations is reached or an
%                          error in the inputs is detected or a NAN or Inf
%                          is detected or a warning is encountered

addpath('C:\Users\wangx\Documents\GitHub\nonlinear-optimization-course\MATLAB\algorithms\modNewton')

% Set dummy values for outputs; prevents errors resulting from bad inputs.  
x        = [];
f        = [];  info.f = f;
f_evals  = 0 ;  info.f_evals = f_evals;
g        = [];  info.g = g;
g_evals  = 0 ;  info.g_evals = g_evals;
H        = [];  info.H = H;
H_evals  = 0 ;  info.H_evals = H_evals;
Hv_evals = 0 ;  info.Hv_evals = Hv_evals;
iter     = 0 ;  info.iter = iter;
status   = 0 ;  info.status = status;

% Check whether enough input parameters were supplied 
if nargin < 2
    fprintf('\n steepest_descent(ERROR):not enough input parameters.\n');
    info.status = -1;
    return
end

% If three input arguments, then define default values for fourth input.
if nargin == 2
    params.dir_type   = 'SteepestDescent';
    params.maxiter    = 1e4;
    params.printlevel = 1;
    params.tol        = 1e-8;
    params.outfileID  = 1;
    params.probname   = '';
end

% Check that the initial point has the desired length.
if length(x0) ~= fun_hands.n
    str = 'x0';
    fprintf('\n backtracking-Armijo(ERROR):Invalid argument %s.\n',str);
    info.status = 1;
    return
end

% Check that all fields required in input params are supplied.
if isfield(params,'dir_type')
    dir_type = params.dir_type;
    if strcmp(dir_type,'SteepestDescent') || strcmp(dir_type,'ModifiedNewton')
    else
        str = 'DIR_TYPE';
        fprintf('\n backtracking-Armijo(ERROR):Invalid argument %s.\n',str);
        info.status = 1;
        return
    end
else
    str = 'DIR_TYPE';
    fprintf('\n backtracking-Armijo(ERROR):control parameter %s not supplied.\n',str);
    info.status = 1;
    return
end

if isfield(params,'maxiter')
    maxiter = params.maxiter;
    if maxiter < 0
        str = 'MAXITER';
        fprintf('\n backtracking-Armijo(ERROR):Invalid control parameter %s.\n',str);
        info.status = 1;
        return
    end
else
    str = 'MAXITER';
    fprintf('\n steepest_descent(ERROR):control parameter %s not supplied.\n',str);
    info.status = 1;
    return
end

if isfield(params,'printlevel')
    printlevel = params.printlevel;
    if printlevel < 0
        str = 'printlevel';
        fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
        info.status = 1;
        return
    end
else
    str = 'PRINTLEVEL';
    fprintf('\n steepest_descent(ERROR):control parameter %s not supplied.\n',str);
    info.status = 1;
    return
end

if isfield(params,'tol')
    tol = params.tol;
    if tol < 0
        str = 'TOL';
        fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
        info.status = 1;
        return
    end
else
    str = 'TOL';
    fprintf('\n steepest_descent(ERROR):control parameter %s not supplied.\n',str);
    info.status = 1;
    return
end
    
if isfield(params,'outfileID')
    outfileID = params.outfileID;
    if outfileID <= 0 || outfileID == 2
        str = 'outfileID';
        fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
        info.status = 1;
        return
    end
else
    outfileID = 1; % standard output (the screen)
end
outfileNAME = fopen(outfileID);

if isfield(params,'probname')
    probname = params.probname;
    if ~ischar(probname)
        str = 'probname';
        fprintf('\n steepest_descen(ERROR):Invalid control parameter %s.\n',str);
        info.status = 1;
        return
    end
else
    probname = ''; 
end

% Check that input fun_hands are all function handles.
f_hand  = fun_hands.f;
g_hand  = fun_hands.g;
H_hand  = fun_hands.H;
Hv_hand = fun_hands.Hv;

if ~isa(f_hand,'function_handle')
    str = 'f_hand';
    fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
    info.status = 1;
    return
end

if ~isa(g_hand,'function_handle')
    str = 'g_hand';
    fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
    info.status = 1;
    return
end

if ~isa(H_hand,'function_handle')
    str = 'H_hand';
    fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
    info.status = 1;
    return
end

if ~isa(Hv_hand,'function_handle')
    str = 'Hv_hand';
    fprintf('\n steepest_descent(ERROR):Invalid control parameter %s.\n',str);
    info.status = 1;
    return
end

% Constant
TINY = eps^(2/3); % Determines if step is too small to make progress.

% Initialization.
x       = x0;
norm_x  = norm(x);
f       = feval(f_hand,x);
g       = feval(g_hand,x);
norm_g  = norm(g);
norm_g0 = norm_g; % Save value at x0 to use in relative stopping condition.
H       = feval(H_hand,x);

% Backtracking Armijo parameters
eta     = 1e-2;
tau     = 0.8;
alpha   = 1;

% Modified Newtons parameters
lam_min = 1e-2;
lam_max = 1e3;

% Save things to make printing easier.
dashedline = repelem('-',1,81) ;
header     = '  Iter         f          norm_g         norm_x        norm_p  ';

% Print column header and value of f at initial point.
if printlevel ~= 0  
  fprintf(outfileID,'%s\n',dashedline);
  fprintf(outfileID,'                            Backtracking Armijo Method \n');
  fprintf(outfileID,'%s\n',dashedline);
  fprintf(outfileID,' maximum iterations    : %g\n',maxiter);
  fprintf(outfileID,' print level           : %g\n',printlevel);
  fprintf(outfileID,' termination tolerance : %1.2e\n',tol);
  fprintf(outfileID,' file for output       : %s\n',outfileNAME);
  fprintf(outfileID,' problem name          : %s\n',probname);
  fprintf(outfileID,' direction method      : %s\n',dir_type);
  fprintf(outfileID,'%s\n',dashedline);
  fprintf(outfileID,'%s\n',header);
  fprintf(outfileID,' %5g %14.7e %14.7e %14.7e', iter, f, norm_g, norm_x);
end


% Main loop: perform the backtracking-Armijo linesearch method.
while (1)
    
    % Check for termination
    if ( norm_g <= tol*max(1,norm_g0) )
       status = 0;
       outcome = ' Relative stopping tolerance reached';
       break
    elseif ( iter >= maxiter )
       status = 1;
       outcome = ' Maximum allowed iterations reached';
       break
    end
    
    % Choose descent direction (only two available)
    if strcmp(dir_type,'SteepestDescent') 
        p = - g;
    else
        [B,~] = modNewton(H, lam_min, lam_max);
        p = - B\g;                 % Solve the Modified Newton direction
    end
    
    % Choose the stepsize taken along the direction p.
    while (1)
        f_new = feval(f_hand,x+alpha*p);
        if f_new <= f+eta*alpha*g'*p
            break
        end
        alpha = alpha*tau;
    end
    
    % Step to actually take.
    norm_p = alpha*norm(p);
    
    % Save norm of current iterate before updating iterate.
    norm_xprev = norm_x;

    % Update the iterate and its associated values.
    iter    = iter + 1;
    x       = x + alpha*p;
    norm_x  = norm(x);
    f       = feval(f_hand,x);
    g       = feval(g_hand,x);
    norm_g  = norm(g);
    H       = feval(H_hand,x);
    f_evals = f_evals + 1;
    g_evals = g_evals + 1;
    H_evals = H_evals + 1;
    
    % Check for NaN/Inf in f, g and H.
    if isnan(f) || isinf(f)
        status = 1;
        outcome = ' ERROR (NaN/Inf when evaluating f)';
        break
    elseif sum(isnan(g(:))) >= 1 || sum(isinf(g(:))) >= 1
        status = 1;
        outcome = ' ERROR (NaN/Inf when evaluating g)';
        break
    elseif sum(isnan(H(:))) >= 1 || sum(isinf(H(:))) >= 1
        status = 1;
        outcome = ' ERROR (NaN/Inf when evaluating H)';
        break
    end

    % Print iterate information, if needed.
    if params.printlevel ~= 0
        fprintf(outfileID,' %14.7e \n', norm_p); % End previous line.
        if mod(iter,20) == 0
          fprintf(outfileID,'%s\n',header); % Print header every 20 iterations.
        end
        fprintf(outfileID,' %5g %14.7e %14.7e %14.7e', iter, f, norm_g, norm_x); % Start next line.
    end

    if norm_p < (1 + norm_xprev)*TINY
        status = 1;
        outcome = ' Step is too small to make additional progress.';
        break
    end
end

% Print termination message, if requested.
if printlevel
  fprintf(outfileID,'\n\n Result      :%s \n', outcome);
  fprintf(outfileID,' Iterations   : %-5g\n', iter);
  fprintf(outfileID,' Final f      : %13.7e\n', f );
  fprintf(outfileID,' Final ||g||  : %13.7e\n', norm_g );
  fprintf(outfileID,' Final ||g0|| : %13.7e\n', norm_g0 );  
  fprintf(outfileID,' Final ||p||  : %13.7e\n', norm_p );  
  fprintf(outfileID,' f_evals      : %-5g\n', f_evals );  
  fprintf(outfileID,' g_evals      : %-5g\n', g_evals ); 
  fprintf(outfileID,' H_evals      : %-5g\n', H_evals ); 
  fprintf(outfileID,'%s\n',dashedline);
end

% Fill output variable inform.
info.f       = f;
info.f_evals = f_evals;
info.g       = g;
info.g_evals = g_evals;  
info.norm_g  = norm_g;
info.H       = H;
info.H_evals = H_evals;
info.iter    = iter;
info.status  = status;

end
function [B, flag] = modNewton(H, lam_min, lam_max)
% The function call aims to compute matrix B for modified Newtons Method
% The arguments are defined as follows:
%
%  Input arguments:
%  ----------------
%    H          : A symmetric matrix.
%    lam_min    : The minimum value of lambda.
%    lam_max    : The maximum limit of lambda.
%
%  Output arguments:
%  -----------------
%    B      : final matrix computed
%    flag   : integer indicating the reason for termination
%               0  Successful termination and no change in the matrix
%               1  Successful termination and change in the matrix
%                      

% Check whether enough input parameters were supplied 
if nargin < 3
    error('\n modNewton(ERROR):not enough input parameters. \n')
end
% Check whether H is symmetric 
if ~isequal(H,H')
    error('\n H should be a symmetric matrix. \n')
end
% Check whether lam_min is greater than zero
if lam_min <= 0
    error('\n lam_min should be greter than zero. \n')
end
% Check whether lam_max is greater than lam_min
if lam_max <= lam_min
    error('\n lam_max should be greater than zero. \n')
end
% Check whether lam_max is beyond the precision limit
if lam_max > 1/eps(2/3)
    error('\n lam_max should be less than the machine precision \n')
end

% Compute Spectral Decomposition
[V, lam] = eig(full(H), 'vector');
mod_lam = zeros(length(lam),1);

% Modify the eigen values
for i = 1:length(lam)
    if abs(lam(i)) > lam_max
        mod_lam(i) = lam_max;
    elseif abs(lam(i)) < lam_min
        mod_lam(i) = lam_min;
    elseif abs(lam(i)) >= lam_min && abs(lam(i)) <= lam_max
        mod_lam(i) = abs(lam(i));
    else
        mod_lam(i) = lam(i);
    end
end

% Obtain the B matrix
B = V*diag(mod_lam)*V';

% Getting the flag
if isequal(lam, mod_lam)
    flag = 0;
else 
    flag = 1;
end

end


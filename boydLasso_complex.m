function [x,z,u,k] = boydLasso_complex(A, b, lambda, rho, alpha,x,u)
% lasso  Solve lasso problem via ADMM
%
% [z, history] = lasso(A, b, lambda, rho, alpha);
%
% Solves the following problem via ADMM:
%
%   minimize 1/2*|| Ax - b ||_2^2 + \lambda || x ||_1
%
% The solution is returned in the vector x.
%
% history is a structure that contains the objective value, the primal and
% dual residual norms, and the tolerances for the primal and dual residual
% norms at each iteration.
%
% rho is the augmented Lagrangian parameter.
%
% alpha is the over-relaxation parameter (typical values for alpha are
% between 1.0 and 1.8).
%
%
% More information can be found in the paper linked at:
% http://www.stanford.edu/~boyd/papers/distr_opt_stat_learning_admm.html
%

t_start = tic;
%Global constants and defaults

MAX_ITER = 1e5;
ABSTOL   = 1e-4;
RELTOL   = 1e-2;
%Data preprocessing

[m, n] = size(A);

% save a matrix-vector multiply
Atb = A'*b;
%ADMM solver

if nargin<6
    x = zeros(n,1);
    z = zeros(n,1);
    u = zeros(n,1);
else
    z=x;
end

% cache the factorization
[L U] = factor(A, rho);


for k = 1:MAX_ITER
    
    % x-update
    q = Atb + rho*(z - u);    % temporary value
    if( m >= n )    % if skinny
        x = U \ (L \ q);
    else            % if fat
        x = q/rho - (A'*(U \ ( L \ (A*q) )))/rho^2;
    end
    
    % z-update with relaxation
    zold = z;
    x_hat = alpha*x + (1 - alpha)*zold;
    z = shrinkage(x_hat + u, lambda/rho);
    
    % u-update
    u = u + (x_hat - z);
    
    % diagnostics, reporting, termination checks
    r_norm = norm(x - z);
    s_norm  = norm(-rho*(z - zold));
    
    eps_pri = sqrt(n)*ABSTOL + RELTOL*max(norm(x), norm(-z));
    eps_dual= sqrt(n)*ABSTOL + RELTOL*norm(rho*u);
    
    
    if (r_norm < eps_pri && s_norm < eps_dual)
        break;
    end
    
end

end

function z = shrinkage(x, kappa)
z = max( 0, 1 - kappa./abs(x) ).*x;
end

function [L U] = factor(A, rho)
[m, n] = size(A);
if ( m >= n )    % if skinny
    L = chol( A'*A + rho*speye(n), 'lower' );
else            % if fat
    L = chol( speye(m) + 1/rho*(A*A'), 'lower' );
end

% force matlab to recognize the upper / lower triangular structure
L = sparse(L);
U = sparse(L');
end
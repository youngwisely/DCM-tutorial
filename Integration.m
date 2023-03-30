function [y,neural]= Integration(P,M,U)
% integrates a MIMO nonlinear system using the Jacobian
% FORMAT [y] = spm_int_J(P,M,U)
% P  - model parameters
% M  - model structure
% U  - input structure or matrix
%
% y  - (v x l)  response y = g(x,u,P)
%__________________________________________________________________________
% Integrates the MIMO system described by
%
%        dx/dt = f(x,u,P,M)
%        y     = g(x,u,P,M)
% or
%        dx/dt = f(x,u,P)
%        y     = g(x,u,P)
%
% using the update scheme:
%
%    x(t + dt) = x(t) + U*dx(t)/dt
%
%            U = (expm(dt*J) - I)*inv(J)
%            J = df/dx
%
% at input times.  This integration scheme evaluates the update matrix (Q)
% at each time point
%
%__________________________________________________________________________
% Copyright (C) 2008 Wellcome Trust Centre for Neuroimaging

% Karl Friston
% $Id: spm_int_J.m 6801 2016-05-29 19:18:06Z karl $


% convert U to U.u if necessary and M(1) to M
%--------------------------------------------------------------------------
dt = U.dt;

% state equation; add [0] states if not specified
%--------------------------------------------------------------------------
f   = fcnchk(M.f,'x','u','P','M');

% and output nonlinearity
%--------------------------------------------------------------------------
g   = @(x,v,P,M)x;

% Initial states and inputs
%--------------------------------------------------------------------------
x   = M.x;

% default delay operator
%--------------------------------------------------------------------------
D = 1;

% integrate
%==========================================================================
for i = 1:size(U.u,1)
    
    % input
    %----------------------------------------------------------------------
        u = U.u(i,:);

    
    % dx(t)/dt and Jacobian df/dx
    %----------------------------------------------------------------------
    if nargout(f) >= 3
        [fx,dfdx,D] = f(x,u,P,M);
        
    elseif nargout(f) == 2
        [fx,dfdx]   = f(x,u,P,M);
        
    else
        fx          = f(x,u,P,M);
        dfdx        = spm_cat(spm_diff(f,x,u,P,M,1));
    end
    
    % update dx = (expm(dt*J) - I)*inv(J)*fx
    %----------------------------------------------------------------------
    x      = spm_unvec(spm_vec(x) + spm_dx(D*dfdx,D*fx,dt),x);
    neural(:,i) = x(:,1);
    
    % output - implement g(x)
    %----------------------------------------------------------------------
    if nargin(g) > 3
        y(:,i) = spm_vec(g(x,u,P,M));
    else
        y(:,i) = spm_vec(g(x,u,P));
    end
    
end

% transpose
%--------------------------------------------------------------------------
y      = real(y');


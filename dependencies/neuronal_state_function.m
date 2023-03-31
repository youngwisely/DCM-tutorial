function [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% State equation for a dynamic [bilinear/nonlinear/Balloon] model of fMRI
% responses
% FORMAT [f,dfdx,D,dfdu] = spm_fx_fmri(x,u,P,M)
% x      - state vector
%   x(:,1) - excitatory neuronal activity            ue
%   x(:,2) - vascular signal                          s
%   x(:,3) - rCBF                                  ln(f)
%   x(:,4) - venous volume                         ln(v)
%   x(:,5) - deoyxHb                               ln(q)
%  [x(:,6) - inhibitory neuronal activity             ui
%
% f      - dx/dt
% dfdx   - df/dx
% dfdu   - df/du
% D      - delays
%
%__________________________________________________________________________
%
% References for hemodynamic & neuronal state equations:
% 1. Buxton RB, Wong EC & Frank LR. t. MRM 39:855-864,
%    1998.
% 2. Friston KJ, Mechelli A, Turner R, Price CJ. Nonlinear responses in
%    fMRI: the Balloon model, Volterra kernels, and other hemodynamics.
%    Neuroimage 12:466-477, 2000.
% 3. Stephan KE, Kasper L, Harrison LM, Daunizeau J, den Ouden HE,
%    Breakspear M, Friston KJ. Nonlinear dynamic causal models for fMRI.
%    Neuroimage 42:649-662, 2008.
% 4. Marreiros AC, Kiebel SJ, Friston KJ. Dynamic causal modelling for
%    fMRI: a two-state model.
%    Neuroimage. 2008 Jan 1;39(1):269-78.
%__________________________________________________________________________
% Copyright (C) 2002-2014 Wellcome Trust Centre for Neuroimaging

% Karl Friston & Klaas Enno Stephan
% $Id: spm_fx_fmri.m 7525 2019-02-05 16:47:50Z peter $


% Neuronal motion
%==========================================================================
P.A   = full(P.A);                       %    linear parameters
P.B   = full(P.B);                       % bi-linear parameters
P.C   = P.C/16;                          % exogenous parameters
P.D   = full(P.D);                       % nonlinear parameters

% implement differential state equation y = dx/dt (neuronal)
%--------------------------------------------------------------------------
f    = x;
  
% one neuronal state per region: diag(A) is a log self-inhibition
%------------------------------------------------------------------
SE     = diag(P.A);
EE     = P.A - diag(exp(SE) + SE);


% flow
%----------------------------------------------------------------------
f(:,1) = EE*x(:,1) + P.C*u(:);
    

% Hemodynamic motion
%==========================================================================

% hemodynamic parameters
%--------------------------------------------------------------------------
%   H(1) - signal decay                                   d(ds/dt)/ds)
%   H(2) - autoregulation                                 d(ds/dt)/df)
%   H(3) - transit time                                   (t0)
%   H(4) - exponent for Fout(v)                           (alpha)
%   H(5) - resting oxygen extraction                      (E0)
%   H(6) - ratio of intra- to extra-vascular components   (epsilon)
%          of the gradient echo signal
%--------------------------------------------------------------------------
H        = [0.64 0.32 2.00 0.32 0.32];

% exponentiation of hemodynamic state variables
%--------------------------------------------------------------------------
x(:,3:5) = exp(x(:,3:5));

% signal decay
%--------------------------------------------------------------------------
sd       = H(1)*exp(P.decay);

% transit time
%--------------------------------------------------------------------------
tt       = H(3)*exp(P.transit);

% Fout = f(v) - outflow
%--------------------------------------------------------------------------
fv       = x(:,4).^(1/H(4));

% e = f(f) - oxygen extraction
%--------------------------------------------------------------------------
ff       = (1 - (1 - H(5)).^(1./x(:,3)))/H(5);


% implement differential state equation f = dx/dt (hemodyna mic)
%--------------------------------------------------------------------------
f(:,2)   = x(:,1) - sd.*x(:,2) - H(2)*(x(:,3) - 1);
f(:,3)   = x(:,2)./x(:,3);
f(:,4)   = (x(:,3) - fv)./(tt.*x(:,4));
f(:,5)   = (ff.*x(:,3) - fv.*x(:,5)./x(:,4))./(tt.*x(:,5));
f        = f(:);


% Neuronal Jacobian
%==========================================================================
[n,m] = size(x);
% one neuronal state per region
%----------------------------------------------------------------------
dfdx{1,1} = EE;
    


% input
%==========================================================================
dfdu{1,1} = P.C;
dfdu{2,1} = sparse(n*(m - 1),length(u(:)));


% Hemodynamic Jacobian
%==========================================================================
dfdx{2,1} = speye(n,n);
dfdx{2,2} = speye(n,n)*(-sd);
dfdx{2,3} = diag(-H(2)*x(:,3));
dfdx{3,2} = diag( 1./x(:,3));
dfdx{3,3} = diag(-x(:,2)./x(:,3));
dfdx{4,3} = diag( x(:,3)./(tt.*x(:,4)));
dfdx{4,4} = diag(-x(:,4).^(1/H(4) - 1)./(tt*H(4)) - (1./x(:,4).*(x(:,3) - x(:,4).^(1/H(4))))./tt);
dfdx{5,3} = diag((x(:,3) + log(1 - H(5)).*(1 - H(5)).^(1./x(:,3)) - x(:,3).*(1 - H(5)).^(1./x(:,3)))./(tt.*x(:,5)*H(5)));
dfdx{5,4} = diag((x(:,4).^(1/H(4) - 1)*(H(4) - 1))./(tt*H(4)));
dfdx{5,5} = diag((x(:,3)./x(:,5)).*((1 - H(5)).^(1./x(:,3)) - 1)./(tt*H(5)));


% concatenate
%--------------------------------------------------------------------------
dfdx      = spm_cat(dfdx);
dfdu      = spm_cat(dfdu);
D         = 1;
function Abunds = SUNSAL(X,EM)
%
% Estimate endmember abundaces using the sparse unmixing via split
% augmented Lagrangian (SUNSAL) method of [1].
% This version of the algorthm is based on the FCLS variant only.
%
%
% Input:
%       X - nData x nVar matrix of observations
%       EM - k x nVar matrix of end member signatures
%
% Output:
%       Abunds - nData x k matrix of the end member abundances
%
%
% Note:
% The is a modification of the SUNSAL funtion written by Jose Bioucas-Dias.
% The original code was downloaded from (http://www.lx.it.pt/~bioucas/code.htm).
% This script has a separate licence agreement.
% All modifications were made by Greig A. Paterson.
%
%
% References:
%
% [1] Bioucas-Dias, J.M., and M. Figueiredo, Alternating direction algorithms
%     for constrained sparse regression: Application to hyperspectral unmixing,
%     in 2nd  IEEE GRSS Workshop on Hyperspectral Image and Signal
%     Processing-WHISPERS' 2010, Raykjavik, Iceland, 2010.
%

%% Do some basic check and input processing

if nargin < 2
    error('SUNSAL:Input', 'At least 2 inputs are required.');
end

if isempty(X)
    error('SUNSAL:Data', 'No data provided.');
end

if isempty(EM)
    error('SUNSAL:EMs', 'No EM data provided.');
end


[nData, nVar] = size(X);
[k, nVar2] = size(EM);

if nVar ~= nVar2
    error('SUNSAL:nVar', 'The number of variables in the data and end member matrices do not match.');
end

% Transponse X to be nVar x nData (to fit the convention of the algorithm)
X = X'; % [nVar x nData]
EM = EM'; % [nVar x k]


%% Set up some default parameters

% maximum number of AL iteration
AL_iters = 500;

% regularization parameter
lambda = 0;
lambda = lambda*ones(k,nData);

% tolerance for the primal and dual residues
tol = 1e-5;

% compute mean norm
norm_y = sqrt(mean(mean(X.^2)));

% rescale M and Y and lambda
EM = EM/norm_y;
X = X/norm_y;
lambda = lambda/norm_y^2;


%% Constants and initializations

mu_AL = 0.01;
mu = 10*mean(lambda(:)) + mu_AL;

[UF,SF] = svd(EM'*EM);
sF = diag(SF);
IF = UF*diag(1./(sF+mu))*UF';

B = ones(1,k);
a = ones(1,nData);

Aux = (IF*B')/inv(B*IF*B');
x_aux = Aux*a;
IF1 = (IF-Aux*B*IF);

yy = EM'*X;

% Initialization
x = IF*EM'*X;
z = x;

% scaled Lagrange Multipliers
d = 0*z;


%% The main body for AL iterations

tol1 = sqrt(k*nData)*tol;
tol2 = sqrt(k*nData)*tol;
ii = 1;
res_p = inf;
res_d = inf;
mu_changed = 0;

while (ii <= AL_iters) && ((abs (res_p) > tol1) || (abs (res_d) > tol2))
    
    % save z to be used later
    if mod(ii,10) == 1
        z0 = z;
    end
    
    % minimize with respect to z
    z =  soft(x-d,lambda/mu);
    
    maskz = (z >= 0);
    z = z.*maskz;
    
    x = IF1*(yy+mu*(z+d))+x_aux;
    
    % Lagrange multipliers update
    d = d -(x-z);
    
    % update mu so to keep primal and dual residuals whithin a factor of 10
    if mod(ii,10) == 1
        
        % primal residue
        res_p = norm(x-z,'fro');
        
        % dual residue
        res_d = mu*norm(z-z0,'fro');
        
        % update mu
        if res_p > 10*res_d
            mu = mu*2;
            d = d/2;
            mu_changed = 1;
        elseif res_d > 10*res_p
            mu = mu/2;
            d = d*2;
            mu_changed = 1;
        end
        
        if  mu_changed
            % update IF and IF1
            IF = UF*diag(1./(sF+mu))*UF';
            Aux = (IF*B')/(B*IF*B');
            x_aux = Aux*a;
            IF1 = (IF-Aux*B*IF);
            mu_changed = 0;
        end
        
    end
    
    ii=ii+1;
    
end

%% Stuff to return

% Transpose the final output to be [nData x k]
Abunds = z';


%% Required subfunctions

function y = soft(x,T)
%
% soft-thresholding function proximity operator for l1 norm
%

if sum(abs(T(:)))==0
    y = x;
else
    y = max(abs(x) - T, 0);
    y = y./(y+T) .* x;
end


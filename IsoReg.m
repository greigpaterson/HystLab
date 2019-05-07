function [yIso] = IsoReg(y)
%
% Function to perform isotonic regression following [1] and [2]. The function
% ensures monotonic increasing data, so the inputs should be appropriately
% pre-processed
%
% Input:
%       y - [n x 1] vector of y values (for magnetics this algorithm
%           assumes that the x values are strictly monotonic)
%
% Output:
%       yIso - [n x 1] of monotonic fitted y values
%
% References:
%      [1] Kruskal, J.B. (1964) Nonmetric multidimensional scaling: a
%          numerical method, Psychometrika 29, 115-129, doi: doi:10.1007/BF02289694.
%
%      [2] Leeuw, J. de, K. Hornik, P. Mair (2009) Isotone Optimization in
%          R: Pool-Adjacent-Violators Algorithm (PAVA) and Active Set
%          Methods, J. Stat. Softw. 32, 1?24. doi: 10.18637/jss.v032.i05.
%
% Last Modified 2019/05/07
%

%% Some input checks

if nargin < 1
    error('IsoReg:Input', 'At least 1 input arguments are required.');
end

if ~isvector(y)
    error('IsoReg:Input', 'y should be a vector.');
end

%% The main function

npts = length(y);
Blocks = 1:npts;

% Intialize monotonic y-values
yIso = y;

ii = 0; % iteration counter (for debugging)

% Find the differences
dY = diff(yIso);
nBad = sum(dY < 0);

% The number of data per block
% This is can be used as a weighting function specified by the user
n_per_block = ones(npts, 1);

while nBad ~= 0
    
    ii = ii + 1;
    
    % Get indices of data in each block
    Inds = cumsum([1; dY>0]);
    Blocks = Inds(Blocks);
    
    % Get the weighted sum of each block
    sumyhat = accumarray(Inds, n_per_block .* yIso);
    
    % Update the number of data per block
    n_per_block = accumarray(Inds, n_per_block);
    
    % Get the mean
    yIso = sumyhat ./ n_per_block;
    
    % Find the differences and the number of non-montoic data
    dY = diff(yIso);
    nBad = sum(dY < 0);
    
end


% Rebuild the monotonic y values
yIso = yIso(Blocks);



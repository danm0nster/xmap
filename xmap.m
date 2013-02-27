function [ X_MY, Y_MX, X1, Y1] = xmap( X, Y, MX, MY, E, tau, L, sampling )
%XMAP calculates cross-mapped estimates Y_MX of Y from MX and vice versa.
%   A library of length L is extracted from the shadow manifolds MX and MY,
%   constructed from X and Y by phase-space reconstruction using E
%   dimensions and a time delay tau. In addition to the cross-mapped
%   estimates X_MY and Y_MX, the original time series points corresponding
%   to the cross-mapped estimates are returned as X1 and Y1.
%
%   XMAP(X, Y, MX, MY, E, tau) computes cross-mapped estimates of two time
%   series X and Y from the phase-space reconstructed shadow manifolds of
%   the other time series (MY and MX), where E is the embedding dimension
%   and tau is the time delay used in the phase-space reconstruction.
%
%   XMAP(X,Y, MX, MY, E, tau, L) computes cross-mapped estimates using a
%   library of L points from the shadow manifolds. If L is not provided, or
%   is empty ([]) then the maximum possible library lentgh is chosen.
%
%   XMAP(X,Y, MX, MY, E, tau, L, sampling) computes cross-mapped estimates,
%   where the parameter sampling is a string whose value determines how the
%   L points are sampled from the shadow manifolds. Possible values are:
%
%     'linear'    The first L points in MX and MY are chosen
%     'random'    A uniformly random set of L points from MX & MY are chosen
%
%   If the parameter sampling is not provided, linear sampling is used.
%
%   Note that the cross-mapped estimates X_MY and Y_MX are time-shifted
%   relative to the original data vectors X and Y. To compare X_MY to X it
%   is necessary to apply the shift X(1+(E-1)*tau:L+(E-1)*tau) to compare
%   with X_MY(1:L).In the case of random sampling the cross-mapped
%   estimates will be shuffled relative to the original input vectors, so
%   time-shifted (for linear sampling) or shuffled (for random sampling)
%   data vectors are returned in the variables X1 and Y1. These can be
%   directly compared to X_MY and Y_MX element by element.
%   
%   see also psembed


%
% Check the input variables and issue an error message if something is
% wrong.
%
MinArg = 6; % Minimum number of required arguments
MaxArg = 8; % Maximum number of allowed arguments
if nargin >= MinArg && nargin <= MaxArg
    if ~(numel(X) == numel(Y))
        error('xmap:argchk','X and Y must have the same length');
    end
    if ~(numel(MX) == numel(MY))
        error('xmap:argchk','MX and MY must have the same size');       
    end    
    if (numel(size(MX)) ~= 2) || (numel(size(MY)) ~= 2)
        error('xmap:argchk','MX and MY must be 2-dimensional arrays');
    end
    %
    % Check that MX and MY have the expected number of dimensions
    %
    [n, ex] = size(MX);
    [n, ey] = size(MY);
    if (ex ~= E) || (ey ~= E)
        error('xmap:argchk','Embedding dimension inconsistent with array size');
    end
    if nargin > MinArg
        if isempty(L)
            L = n; % Set L to maximum possible value
        elseif L > n
            error('xmap:argchk','L exceeds data size');
        end
    else
        L = n; % Set L to maximum possible value
    end
    if nargin > MinArg+1
        if ~ischar(sampling)
            error('xmap:argchk','sampling must be a string');
        end
    else
        sampling = 'linear'; % Default if no option is given
    end
elseif nargin > MaxArg
    error('xmap:argchk','Too many arguments');
elseif nargin < MinArg
    error('xmap:argchk','Too few arguments');
end

%
% Extract a library of 'L' points from MX and MY
%
switch (sampling)
    case 'linear'
        MX = MX(1:L,:);
        MY = MY(1:L,:);
    case 'random'
        % Uniformly distributed random intergers between 1 and n
        idx = randi(n,size([1:L]')); 
        MX = MX(idx,:);     
        idy = randi(n,size([1:L]'));
        MY = MY(idy,:);        
    otherwise
        errorstring = sprintf('Unknown sampling type: %s',sampling);
        error('xmap:sampling', errorstring);
end
%
% Find the E+1 nearest neighbors of every point in MX and MY. The closest
% neighbor to MX(j,:) is itself, so E+2 neighbors are found.
%
[n,d]=knnsearch(MX,MX,'k',E+2,'distance','euclidean');
%
% n contains the indices to MX for the E+2 nearest nighbors (the first one
% being each point in MX itself). It has E+2 columns, but we only need
% columns 2:E+2, which are the indices of the E+1 neighbors. To convert
% indices to MX to indices to Y (or X) we apply the formula:
%
%  k = n + (E-1)*tau,
%
% where n is an index to MX and k an index to Y. So now the Y-values
% corresponding to the nearest neighbors to points on the shadow manifold
% MX are given by the expression below. Each point in MX is represented as
% a row, and the E+1 columns are the corresponding Y-values needed to
% estimate Y at the time concurrent with the point in MX.
%
switch (sampling)
    case 'linear'
        data_nn = Y(n(:,2:E+2)+(E-1)*tau);        
    case 'random'
        data_nn = Y(idx(n(:,2:E+2))+(E-1)*tau);       
end

%
% Now the weights can be computed as exp( -d(j)/(d(1) ) for j=1:E+1. These
% distances for the point with index p are given by d(p,2:E+2).
% 
% I add the small constant EPS to the denominator in order to avoid
% floating point overflow, resulting in NaN. Perhaps a better choice would
% be to let it depend on the numerator. Investigate and test this further!
%
EPS = 1e-9;
w = zeros(L,E+1);
for p=1:L
    w(p,:) = exp(-d(p,2:E+2)/(EPS+d(p,2)));
    w(p,:) = w(p,:)/sum(w(p,:));
end
%
% Compute cross-mapped estimates as an exponentially weighted sum
%
Y_MX = w .* data_nn;
Y_MX = sum(Y_MX,2);
%
% And now we can do everything in the other direction, i.e. estimating X
% from MY.
%
[n,d]=knnsearch(MY,MY,'k',E+2,'distance','euclidean');
switch (sampling)
    case 'linear'
        data_nn = X(n(:,2:E+2)+(E-1)*tau);        
    case 'random'
        data_nn = X(idy(n(:,2:E+2))+(E-1)*tau);        
end
for p=1:L
    w(p,:) = exp(-d(p,2:E+2)/(EPS+d(p,2)));
    w(p,:) = w(p,:)/sum(w(p,:));
end
X_MY = w .* data_nn;
X_MY = sum(X_MY,2);
%
% Construct comparison vectors from the original data vectors
%
switch (sampling)
    case 'linear'
        X1 = X(1+(E-1)*tau:L+(E-1)*tau);
        Y1 = Y(1+(E-1)*tau:L+(E-1)*tau);
    case 'random'
        X1 = X(idy+(E-1)*tau);
        Y1 = Y(idx+(E-1)*tau);
end
end


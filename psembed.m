function [ M ] = psembed( tsdata, E, tau)
%PSEMBED Phase space embedding of a times series. 
%   PSEMBED(tsdata, E, tau) performs a phase space reconstruction
%   (embedding) in E dimensions, with a time lag tau, based on the time
%   series in tsdata. 
%
%   see also xmap

switch (nargin)
    case 3
        %
    otherwise
        error('psembed:argchk','Wrong number of arguments arguments');
end
%
% N is the data size, i.e. the number of points in tsdata
%
N = numel(tsdata);
%
% N_vec is the number of phase space vectors in the embedding, using the N
% points from tsdata to embed in E dimensions with time lag tau.
%
N_vec = N-(E-1)*tau;
%
% Allocate memory for the embedded manifold in the matrix M.
%
M = zeros(N_vec,E);
%
% Now construct the phase-space embedded manifold, using the time-delayed
% values from tsdata.
%   
for m=1:N_vec
    for j=1:E
        M(m,j) = tsdata(m+(E-j)*tau);
    end
end
end


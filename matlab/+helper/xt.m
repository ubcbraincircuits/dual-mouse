function t = xt(X,fs,DIM)
% Generate time vector
%
% Inputs:
%   X            (required, time series data)
%   fs           (required, sampling frequency)
%   DIM          (optional, dimension to calculate time across. Assume 
%                assumes dimension with highest length is time)
%
% Outputs:
%   t            (array of sample times)
% 
% Usage:  t = xt(X,fs,DIM)

if nargin < 3 || isempty(DIM), [~,DIM] = max(size(X)); end

t = (0:size(X,DIM)-1)/fs;

end
function h = isnormal(data)
% returns logical 
%   true if data is normally distributed

    h = ~kstest(zscore(data));
end
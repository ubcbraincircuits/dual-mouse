function uboxplot(varargin)
% boxplot of unequal sizes. fill in with nans. must be column vectors input

N = 1;
for i = 1:nargin
    N = max([N, size(varargin{i}, 1)]);
end


tmp = nan(N, nargin);
for i = 1:nargin
    tmp(1:size(varargin{i}),i) = varargin{i};
end

boxplot(tmp);
end
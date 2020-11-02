function B = processBehavior(pathname, filenames, gaussParams, stdmult)
if nargin < 3 || isempty(gaussParams), gaussParams = [5, 100]; end
if nargin < 4 || isempty(stdmult), stdmult = 1; end

for i = 1:length(filenames)
    load([pathname, filenames{i}], 'left_whisk', 'right_whisk', ...
        'i1', 'i2', 't1', 't2', ... 
        'left_forelimb', 'right_forelimb')

    B.times(i,:) = [i1, i2, t1, t2];        

    % Get behaviour vectors
    B.whisk.left{i} = getBehIdx(left_whisk, gaussParams, stdmult);
    B.whisk.right{i} = getBehIdx(right_whisk, gaussParams, stdmult);
    B.FL.left{i} = getBehIdx(left_forelimb, gaussParams, stdmult);
    B.FL.right{i} = getBehIdx(right_forelimb, gaussParams, stdmult);


    % Create exclusive behaviours      
    whiskLeft = B.whisk.left{i} - B.FL.left{i};
    whiskLeft = logical(whiskLeft > 0);
    whiskRight = B.whisk.right{i} - B.FL.right{i}; 
    whiskRight = logical(whiskRight > 0);

    flLeft = B.FL.left{i} - B.whisk.left{i};
    flLeft = logical(flLeft > 0);
    flRight = B.FL.right{i} - B.whisk.right{i};
    flRight = logical(flRight > 0);

    % Get trial phase indices
    bIdx = zeros(size(whiskLeft)); bIdx(1:t1) = 1; % before
    dIdx = zeros(size(whiskLeft)); dIdx(i1:i2) = 1; % during
    aIdx = zeros(size(whiskLeft)); aIdx(t2:end) = 1; % after


    B.exclusive.whiskLeftBefore{i} = whiskLeft .* bIdx;
    B.exclusive.whiskLeftDuring{i} = whiskLeft .* dIdx;
    B.exclusive.whiskRightDuring{i} = whiskRight .* dIdx;
    B.exclusive.flLeftBefore{i} = flLeft .* bIdx;
    B.exclusive.flLeftDuring{i} = flLeft .* dIdx;
    B.exclusive.flRightDuring{i} = flRight .* dIdx;
    B.exclusive.whiskLeftAfter{i} = whiskLeft .* aIdx;
end
end


function idx = getBehIdx(Bvec, gaussParams, stdmult, lenBrain)
% Get binary behaviour vector
%   Ensure length of brain data >= length of behaviour data
%   Smooth behaviour gradient with Gaussian kernel
%   Threshold at mean + 1 sd
if nargin < 2 || isempty(gaussParams), gaussParams = [5, 100]; end
if nargin < 3 || isempty(stdmult), stdmult = 1; end
if nargin < 4 || isempty(lenBrain), lenBrain = inf; end

if length(Bvec) > lenBrain
    Bvec = Bvec(1:lenBrain - 1);
end

smB = helper.gauss1(Bvec, gaussParams(1), gaussParams(2));
idx = logical(smB > mean(smB) + (stdmult*std(smB)));
end
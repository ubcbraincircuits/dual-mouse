function [tform, regFrames] = getTforms(data, refIdx, plotFigs, verbose, cpIndex)
% get registration transformations for all frames in a stack
%
% Input:
%   data:       (HxWxN, where [HxW] are image size & N is number of frames)
%   refIdx:     (index of reference frame in stack: refIdx <= N)
%   plotFigs:   (Boolean: True to show registered maps in tiled image)
%   verbose:    (Boolean: True to print progress to command window)
%   cpIndex:    (Index for using control point selection method,
%                otherwise, use default MATLAB registration)
%
% Output:
%   tform:      (Cell array of transformations: length N)
%   regFrames:  (Tiled image of registered maps for visualization)

if nargin < 3 || isempty(refIdx), refIdx = chooseRef(data); end
if nargin < 3 || isempty(plotFigs), plotFigs = 1; end
if nargin < 4 || isempty(verbose), verbose = 1; end
if nargin < 5 || isempty(cpIndex), cpIndex = []; end


N = size(data, 3);
regFrames = zeros(size(data));
tform = cell(1,N);

[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumStepLength = 6.25e-4;

if verbose
    disp(['Reference index: ', num2str(refIdx)])
end

fixed = data(:, :, refIdx);
fixedtmp = imadjust(helper.MM_flat(fixed));
for i = 1:N
    if i ~= refIdx
        moving = data(:, :, i);
%         moving = imhistmatch(moving, fixed);
        movingtmp = imadjust(helper.MM_flat(moving));

        % Use control point registration when needed
        if ismember(i, cpIndex)
            [tform{i}, movingRegistered] = cpregister(moving, fixed);
                regFrames(:, :, i) = movingRegistered;
        else
            tform{i} = imregtform(movingtmp, fixedtmp, 'similarity', ...
                optimizer, metric);
%             tformEstimate = imregtform(movingtmp, fixedtmp, 'rigid', ...
%                 optimizer, metric);
%             
%             tform{i} = imregtform(movingtmp, fixedtmp, 'similarity', ...
%                 optimizer, metric, 'InitialTransformation', tformEstimate);

            regFrames(:, :, i) = imwarp(moving, tform{i}, ...
                        'OutputView', imref2d(size(moving)));

        end
    
    else
        tform{i} = [];
        regFrames(:, :, i) = fixed;
    end
    
    if verbose
        disp([num2str(i),'/',num2str(N),' complete'])
    end
end


if plotFigs
    figure
    imagesc(helper.tileMaps(regFrames))
    colormap gray
    grid on
    title('Registered Frames')
end

end


function refIdx = chooseRef(data)

[H, W, ~] = size(data);

% plot all maps and choose a reference
tiledImage = tileMaps(data);
numCols = size(tiledImage, 2) / W;

f = figure; 
imagesc(tiledImage), colormap gray
[x, y] = ginput(1); 
close(f);

x = ceil(x/W);
y = ceil(y/H);

refIdx = x + (numCols*(y-1));
end


function [tform, movingRegistered] = cpregister(moving, fixed)

[mp,fp] = cpselect(moving,fixed,'Wait',true);
tformEstimate = fitgeotrans(mp,fp,'NonreflectiveSimilarity');


[optimizer, metric] = imregconfig('monomodal');
optimizer.MaximumStepLength = 6.25e-5;
optimizer.MaximumIterations = 100;

tform = imregtform(moving, fixed, 'similarity', optimizer, metric, ...
    'InitialTransformation', tformEstimate);
movingRegistered = imwarp(moving,tform,'OutputView',imref2d(size(fixed)));


end

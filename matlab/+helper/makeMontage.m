function M = makeMontage(image_data, fs, tvec)
% Make a montage of brain maps from image data
% Inputs:
%   image_data of shape Height x Width x Frames
%   fs - sampling rate
%   tvec, time vector specifying timepoints to show image data
%
% Output:
%   2D array of concatenated frames from image_data at timepoints specified
%   in tvec


t = xt(image_data, fs, 3); % time vector for all frames in image_data

% if tvec is not specified, make montage with 10 evenly spaced frames
if nargin < 3 || isempty(tvec)
    tvec = linspace(t(1), t(end), 10);
end

% make the montage
M = [];
for i = 1:length(tvec)
    [~,idx] = min(abs(t-tvec(i)));
    M = [M, image_data(:,:,idx)];
end
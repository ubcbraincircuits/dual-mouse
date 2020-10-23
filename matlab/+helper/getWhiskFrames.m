function wFrames = getWhiskFrames(data, wEvents, window)
% Function to get temporal signals surrounding whisk initiation. Returns a
% 4D array of brain maps vs time for each whisk event.
%
% Input:
%   data        (brain data of shape HxWxT)
%   wEvents     (index of whisk frames)
%   window      (frames before and after whisk initiation)
%
% Output:
%   wFrames (4D array of shape H x W x T x N, where N is number of whisks)


[H, W, ~] = size(data);
wFrames = zeros(H, W, 2*round(window/2) + 1);
win = round(window/2);

numWhiskEvents = numel(wEvents);
for j = 1:numWhiskEvents
    wFrames(:,:,:,j) = data(:, :, wEvents(j)-win: wEvents(j)+win);
end

end
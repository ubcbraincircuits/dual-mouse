function wEvents = getWhiskEvents(idx, window)
% Function to get whisking events. Takes a thresholded behavior gradient
% and returns binary array of whisk initiation events.
%
% Input:
%   idx      (thresholded behavior gradient)
%   window   (frames surrounding whisk initiation event)
% Output:
%   wEvents  (array of whisk initiation indexes)


% create array of whisk initiation events
wEvents = zeros(size(idx));
for i = 1:length(idx)-1
    if idx(i) == 0 && idx(i+1) == 1
        wEvents(i+1) = 1;
    end
end

% remove events that occur within one window length of each other
wEventIdx = find(wEvents == 1);
for i = 1:length(wEventIdx)-1
    
    idx1 = wEventIdx(i);
    idx2 = wEventIdx(i+1);
    
    if wEventIdx(i) - round(window/2) <= 0
        wEvents(idx1) = 0;
    elseif wEventIdx(i+1) - wEventIdx(i) < window
        wEvents(idx2) = 0;
    elseif wEventIdx(i+1) + round(window/2) > length(idx)
        wEvents(idx2) = 0;
    end
end

wEvents = find(wEvents == 1);

end
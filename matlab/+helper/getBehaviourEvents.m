function [bOn, bOff] = getBehaviourEvents(behaviour)
bOn = find(diff(behaviour) == 1);
bOff = find(diff(behaviour) == -1);
% if length(bOn) > length(bOff)
%     % add ending index
%     bOff = [bOff, length(behaviour)];
% 
%     % remove last start index
% %     bOn(end) = [];
% elseif length(bOn) < length(bOff)
% %     % add starting index
%     bOn = [1, bOn];
% 
%     % remove first ending index
% %     bOff(1) = [];
% end


while bOff(1) < bOn(1)
    bOff(1) = [];
end
while bOn(end) > bOff(end)
    bOn(end) = [];
end

bOn = squeeze(bOn);
bOff = squeeze(bOff);
end
function [bOn, bOff] = getBehaviourEvents(behaviour, window)
bOn = find(diff(behaviour) == 1);
bOff = find(diff(behaviour) == -1);

if numel(bOff) > 0 && numel(bOn) > 0

while bOff(1) < bOn(1)
    bOff(1) = [];
end
while bOn(end) > bOff(end)
    bOn(end) = [];
end

bOn = squeeze(bOn);
bOff = squeeze(bOff);
end



% remove events that occur within one window length of each other
bOnIdx = ones(size(bOn));
rw = round(window/2);
for i = 1:length(bOn)-1
    if (bOn(i) - rw <= 0) || (bOn(i+1) - bOn(i) < window) || (bOn(i) + rw > length(behaviour) )
        bOnIdx(i) = 0;
    end
end

bOn = bOn(logical(bOnIdx));
bOff = bOff(logical(bOnIdx));

end
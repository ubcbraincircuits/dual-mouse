function bEnd = limitBLength(behavior_array, brain_video)

d = size(brain_video,3);

if length(behavior_array) > d
    bEnd = d-1;
else
    bEnd = length(behavior_array);
end

end
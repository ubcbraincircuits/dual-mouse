function reg_maps = applyTform(frames, tform)

reg_maps = frames;

for i = 1:length(tform)
    if isempty(tform{i}), continue; end
        reg_maps(:,:,i) = imwarp(frames(:,:,i), tform{i}, ...
            'OutputView', imref2d(size(frames(:,:,i))));
    disp([num2str(i),'/',num2str(length(tform))])
end

end
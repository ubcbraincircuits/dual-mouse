function reg_maps = applyTform(frames, tform)

reg_maps = frames;

for i = 1:length(tform)
    try
        reg_maps(:,:,i) = imwarp(frames(:,:,i), tform{i}, ...
            'OutputView', imref2d(size(frames(:,:,:,i))));
    catch
    end
    disp([num2str(i),'/',num2str(length(tform))])
end

end
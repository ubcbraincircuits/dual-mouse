function frames = getSingleFrames(filePath, fileList, H, W, both)
if nargin < 5 || isempty(both), both = true; end


N = length(fileList);
if both
    frames = uint8(zeros(H, W, N*2));
else
    frames = uint8(zeros(H, W, N));
end

for i = 1:N
    load([filePath, fileList{i}], 'left_frame')
    frames(:,:,i) = imrotate(left_frame, 90);
    
    if both
        load([filePath, fileList{i}], 'right_frame')
        frames(:,:,N+i) = imrotate(right_frame, -90);
    end
end


end
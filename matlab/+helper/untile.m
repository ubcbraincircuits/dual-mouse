function stack = untile(tiledImage, imageSize)
if numel(imageSize) == 2
    H = imageSize(1); W = imageSize(2);
elseif numel(imageSize) == 1
    H = imageSize; W = imageSize;
else
    error('imageSize must be two element array denoting height and width');
end

numRows = size(tiledImage,1) / H;
numCols = size(tiledImage,2) / W;

stack = zeros(H, W, numRows * numCols);
for i = 1:numRows
    for j = 1:numCols
        stack(:,:,j+numCols*(i-1)) = tiledImage( (H*(i-1)+1):(H*i), ...
            (W*(j-1)+1):(W*j) );
    end
end

end
function tiledImage = tileMaps(data, normalization)
% Makes a tiled image of brain maps for viewing
%
% Input: data (HxWxN, where [H,W] are image size & N is number of images)
% Output: tiledImage

if nargin<2, normalization = 'None'; end

factorList = [];
[H, W, N] = size(data);
for i = 1:N
    if rem(N, i) == 0
        factorList = cat(1, factorList, i);
    end
end


numRows = factorList(round(length(factorList)/2));
numCols = N/numRows;

tiledImage = zeros(H*numRows, W*numCols);

for i = 1:numRows
    for j = 1:numCols
        try
            img = data(:,:,j+numCols*(i-1));
            switch normalization
                case 'None'
                otherwise
                    img = (img - min(img(:)))./(max(img(:)) - min(img(:)));
            end
            tiledImage( (H*(i-1)+1):(H*i), ...
                (W*(j-1)+1):(W*j) ) = img;
        catch
            tiledImage( (H*(i-1)+1):(H*i), ...
                (W*(j-1)+1):(W*j) ) = zeros(H, W);            
        end
    end
end

end
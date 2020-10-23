function roi = draw_roi(img,numroi,fillvalue)
% Draw regions of interest on image data. This script was adapted from 
% https://www.mathworks.com/matlabcentral/answers/uploaded_files/113201/draw_multiple_polygons.m
% 
% Inputs:
%   img          (required, image data - dimensions: Height x Width x Time)
%   numroi       (optional, number of ROIs - default = 1)
%
% Outputs:
%   roi          (roi mask)
% 
% Usage:  roi = draw_roi(img, numroi);

if nargin < 2 || isempty(numroi), numroi = 2; end
if nargin < 3 || isempty(fillvalue), fillvalue = nan; end

roihandle = figure();
subplot(1,2,1), imagesc(mean(img,3)), title('Original')
subplot(1,2,2), imagesc(mean(img,3)), title('ROIs overlaid'), 
colormap gray, hold on

% Ask user to draw freehand mask.
message = sprintf(['Left click to draw vertices in the left image.',...
    '\nRight click to finish. Vertices can added by pressing A.',...
    '\nDouble click in the middle to accept.']);

answer = questdlg(message,'','Quit','Continue','Continue');

switch answer
    case 'Continue'
        regionCount = 0;
        roi = false(size(img));
        while regionCount < numroi
            regionCount = regionCount + 1;

            subplot(1, 2, 1); % Switch to img axes.
            [single_roi, xi, yi] = roipoly();

            % Draw the polygon over the img on the right.
            subplot(1, 2, 2);
            plot(xi, yi, 'r-', 'LineWidth', 2);

            % Create the combined roi mask
            roi = roi | single_roi;
        end
        roi = single(roi);
        close(roihandle)
    case 'Quit'
        roi = ones(size(img),'single');
end

roi(roi == 0) = fillvalue;

end
function timeseries = getTimeseries(data, CL, CR, radius)
% Get ROI timeseries from image stack. If no ROIs are specified return 
% global median
%
% Inputs:
%   data        (required, image stack with frames along the 3rd dimension)
%   CL, CR      (optional, ROI coordinates in format from get_coords)
%   radius      (optional, px surrounding coordinate location)


if nargin == 1
    timeseries = squeeze(nanmedian(nanmedian(data,2),1));
else
    
    T = size(data, 3);
    R = size(CL, 1);
    timeseries = zeros(T, R);
    
    for i = 1:R
        regiony  = CL(i,2)-radius:CL(i,2)+radius;
        regionxL = CL(i,1)-radius:CL(i,1)+radius;
        regionxR = CR(i,1)-radius:CR(i,1)+radius;
        timeseries(:, i) = squeeze(nanmedian(nanmedian(...
                    data(regiony,[regionxL,regionxR],:), ...
                    2),1));
    end
    

end

end
function filtimage = MM_flat(imagedata,filterLargeDia,filterSmallDia,scaleimage)
% Flatten image illumination using FFT bandpass filter
% Adapted from ImageJ FFT bandpass filter by Joachim Walter
% original -- https://imagej.nih.gov/ij/plugins/download/FFT_Filter.java
%
%
%
% 

if nargin<2, filterLargeDia = 40; end
if nargin<3, filterSmallDia = 1; end
if nargin<4, scaleimage = 4; end

imagedata = imresize(imagedata,scaleimage);
origsize = size(imagedata);
maxN = max(origsize);


i=2;
while i<1.5*maxN
    i = 2*i;
end

filterLarge = 2.0*filterLargeDia / i;
filterSmall = 2.0*filterSmallDia / i;


ip2 = cat(2,imagedata,fliplr(imagedata),imagedata);
ip2 = cat(1,ip2,flipud(ip2));
ip2 = ip2(1:i,1:i);

% take discrete Hartley transform and apply bandpass filter
fht = real(fft2(ip2)) - imag(fft2(ip2));
fftfilter = filterLargeSmall(fht,filterLarge,filterSmall);
filtimage = ifft2(fftfilter);
filtimage = real(filtimage(1:origsize(1),1:origsize(2)))-imag(filtimage(1:origsize(1),1:origsize(2)));
filtimage = imresize(filtimage,1/scaleimage);

end


% bandpass filter
function fht = filterLargeSmall(ip2, filterLarge, filterSmall)

maxN = max(size(ip2));
fht = ip2;
filt = ones(maxN*maxN,1);

scaleLarge = filterLarge*filterLarge;
scaleSmall = filterSmall*filterSmall;

% loop over rows
for j = 1 : maxN/2
    row = j * maxN;
    backrow = (maxN-j)*maxN;
    rowFactLarge = exp(-(j*j) * scaleLarge);
    rowFactSmall = exp(-(j*j) * scaleSmall);
    
    % loop over columns
    for col = 1 : maxN/2
        backcol = maxN-col;
        colFactLarge = exp(-(col*col) * scaleLarge);
        colFactSmall = exp(-(col*col) * scaleSmall);
        
        factor = (1-rowFactLarge*colFactLarge) * rowFactSmall*colFactSmall;
        
        fht(col+row) = fht(col+row) * factor;
        fht(col+backrow) = fht(col+backrow) * factor;
        fht(backcol+row) = fht(backcol+row) * factor;
        fht(backcol+backrow) = fht(backcol+backrow) * factor;
        
        filt(col+row) = filt(col+row) * factor;
        filt(col+backrow) = filt(col+backrow) * factor;
        filt(backcol+row) = filt(backcol+row) * factor;
        filt(backcol+backrow) = filt(backcol+backrow) * factor;
    end
end

rowmid = maxN * (maxN/2);
rowFactLarge = exp(- (maxN/2)*(maxN/2) * scaleLarge);
rowFactSmall = exp(- (maxN/2)*(maxN/2) * scaleSmall);

fht(maxN/2) = fht(maxN/2) * (1 - rowFactLarge) * rowFactSmall; % (maxN/2,0)
fht(rowmid) = fht(rowmid) * (1 - rowFactLarge) * rowFactSmall; % (0,maxN/2)
fht(maxN/2 + rowmid) = fht(maxN/2 + rowmid) * (1 - rowFactLarge*rowFactLarge) * rowFactSmall*rowFactSmall; % (maxN/2,maxN/2)

filt(maxN/2) = filt(maxN/2) * (1 - rowFactLarge) * rowFactSmall; % (maxN/2,0)
filt(rowmid) = filt(rowmid) * (1 - rowFactLarge) * rowFactSmall; % (0,maxN/2)
filt(maxN/2 + rowmid) = filt(maxN/2 + rowmid) * (1 - rowFactLarge*rowFactLarge) * rowFactSmall*rowFactSmall; % (maxN/2,maxN/2)


% loop along row 0 and maxN/2
rowFactLarge = exp(- (maxN/2)*(maxN/2) * scaleLarge);
rowFactSmall = exp(- (maxN/2)*(maxN/2) * scaleSmall);

for col = 1 : maxN/2
    backcol = maxN-col;
    colFactLarge = exp(- (col*col) * scaleLarge);
    colFactSmall = exp(- (col*col) * scaleSmall);
    
    fht(col) = fht(col) * (1-colFactLarge) * colFactSmall;
    fht(backcol) = fht(backcol) * (1-colFactLarge) * colFactSmall;
    fht(col+rowmid) = fht(col+rowmid) * (1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    fht(backcol+rowmid) = fht(backcol+rowmid) * (1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    filt(col) = filt(col) * (1-colFactLarge) * colFactSmall;
    filt(backcol) = filt(backcol) * (1-colFactLarge) * colFactSmall;
    filt(col+rowmid) = filt(col+rowmid) * (1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    filt(backcol+rowmid) = filt(backcol+rowmid)*(1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;	
end

% loop along column 0 and maxN/2
colFactLarge = exp(- (maxN/2)*(maxN/2) * scaleLarge);
colFactSmall = exp(- (maxN/2)*(maxN/2) * scaleSmall);

for j = 1 : maxN/2
    row = j*maxN;
    backrow = (maxN-j)*maxN;
    rowFactLarge = exp(- (j*j) * scaleLarge);
    rowFactSmall = exp(- (j*j) * scaleSmall);
    
    fht(row) = fht(row)*(1-rowFactLarge)*rowFactSmall;
    fht(backrow) = fht(backrow)*(1-rowFactLarge)*rowFactSmall;
    fht(row+maxN/2) = fht(row+maxN/2)*(1-rowFactLarge*colFactLarge)*rowFactSmall*colFactSmall;
    fht(backrow+maxN/2) = fht(backcol+rowmid) * (1-rowFactLarge*colFactLarge) * colFactSmall*rowFactSmall;
    filt(row) = filt(row) * (1-rowFactLarge) * rowFactSmall;
    filt(backrow) = filt(backrow) * (1-rowFactLarge) * rowFactSmall;
    filt(row+maxN/2) = filt(row+maxN/2) * (1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;
    filt(backrow+maxN/2) = filt(backrow+maxN/2)*(1-colFactLarge*rowFactLarge) * colFactSmall*rowFactSmall;	
    
end


end
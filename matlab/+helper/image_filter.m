function filt = image_filter(I)
% assume time on 3rd dimension


fs = 28.815;


[b,a] = butter(2,[0.01 12]/(fs/2), 'bandpass');
% [b,a] = butter(2,[0.4 4]/(fs/2), 'bandpass');
I = permute(I, [3, 1, 2]);
filt = filtfilt(b,a,I);
filt = permute(filt, [2, 3, 1]);


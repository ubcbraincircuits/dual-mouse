function I = open_raw(filename)
resolution = [256, 256];

file = dir(filename);
fid = fopen([file.folder,'\',file.name], 'r');

if fid == -1, error(['Can not open ', file.name]); end


disp(['Opening ', filename,' ...'])
I = fread(fid, file.bytes, 'single');
fclose(fid);


I = reshape(I,resolution(1),resolution(2),[]);
% I = permute(I,[2 3 4 1]);


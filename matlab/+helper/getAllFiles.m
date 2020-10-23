function fileList = getAllFiles(dirName, ext)
% return list of files from specified directory 
%
% input:
%   dirName (directory name, string)
%   ext     (filename extension, string, eg .tif, .mat, .csv, etc)
%
% output:
%   fileList (list of file names, cell array of strings)

if nargin < 2 || isempty(ext), ext = '.'; end

assert(contains(ext, '.'), ...
    'Filename extension must include the dot (eg .mat or .csv).');

dirData = dir(dirName); 

% exclude index for directories
dirIndex = [dirData.isdir]; 
fileList = {dirData(~dirIndex).name}'; 

% restrict fileList to specified filename extension
fileList = fileList( contains(fileList, ext) );

end
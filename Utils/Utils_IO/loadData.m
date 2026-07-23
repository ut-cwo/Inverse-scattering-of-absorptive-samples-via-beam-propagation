function [data, params] = loadData(filePath, fileName)
%% loadFieldData: Load either a .mat or .h5 file.
% Arguments:
%   filePath - char Array to directory of data
%   fileName - char Array with sample name
%
% Returns:
%   data - struct containing loaded contents
%   params - data params
arguments
    filePath    (1,:) char
    fileName    (1,:) char
end

% Define full path
fullPath = [filePath fileName];

if ~isfile(fullPath)
    error('File not found: %s', fullPath);
end

[~,~,ext] = fileparts(fullPath);
ext = lower(ext);

switch ext
    case '.mat'
        % Standard MATLAB file
        data = load(fullPath);
    case {'.h5', '.hdf5'}
        % H5/HDF5 loading
        [data, params]   = largeFileIO.readAllChunks(filePath, fileName, 1);
        data        = single(data)/255; 
    otherwise
        error('Unsupported file type: %s', ext);
end
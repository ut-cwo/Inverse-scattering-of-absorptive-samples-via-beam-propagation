%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                              loadData.m                                 %
% _______________________________________________________________________ %
%                     Data loading helper function                        %
%                                                                         %
%                                                                         %
%                                         July 22, 2026 by Peter Wagenaar %
%                                                                         %
% References:                                                             %
% [1]   Wagenaar, Peter, et al. "Inverse-scattering of absorptive samples %
%       via beam propagation." bioRxiv (2026).                            %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/ut-cwo/Inverse-scattering-in-biological-   %
%               samples-via-beam-propagation                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function data = loadData(filePath, fileName)
%% loadFieldData: Load a .mat file.
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
        data            = load(fullPath);
    otherwise
        error('Unsupported file type: %s', ext);
end
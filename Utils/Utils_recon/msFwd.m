%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                msFwd.m                                  %
% _______________________________________________________________________ %
%              Forward model helper function for mainMSBP.m               %
%                                                                         %
% This function is a reformulation of Multislice_Forward from previous    %
% versions of the code. It has been reworked to increase execution speed  %
% and improve memory utilization.                                         %
%                                                                         %
%                                         July 22, 2026 by Peter Wagenaar %
%                                          August 16, 2025 by Jeongsu Kim %
%                                   July 25, 2020 by Shwetadwip Chowdhury %
%                                                                         %
% References:                                                             %
% [1]   Wagenaar, Peter, et al. "Inverse-scattering of absorptive samples %
%       via beam propagation." bioRxiv (2026).                            %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/ut-cwo/Inverse-scattering-in-biological-   %
%               samples-via-beam-propagation                              %
%                                                                         %
%                                                                         %
% [2]   Kim, Jeongsoo, et al. "Inverse-scattering in biological samples   %
%       via beam-propagation." bioRxiv (2025).                            %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/ut-cwo/Inverse-scattering-of-absorptive-   %
%               samples-via-beam-propagation                              %
%                                                                         %
%                                                                         %
% [3]   S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and   %
%       L. Waller, "High-resolution 3D refractive index microscopy of     %
%       multiple-scattering samples from intensity images,"               %
%       Optica 6, 1211-1219 (2019)                                        %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/Waller-Lab/multi-slice                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [arr] = msFwd(p, r, obj, fU_in, arr)
%% msFwd: Simulate field measurement using multilsice beam propagation
% Arguments:
%   p         - System Parameters Struct
%   r         - Reconstruction Parameters Struct
%   obj       - 3D RI volume
%   fU_in     - 2D field of incoming plane wave
%   arr       - Struct holding pre-allocated field volume
% Returns:
%   arr       - Struct storing 2D/3D field allocations
arguments
    p           struct
    r           struct
    obj         (:,:,:) {mustBeUnderlyingType(obj,'single')}
    fU_in       (:,:)   {mustBeUnderlyingType(fU_in,'single')}
    arr         struct
end

% Incident Plane Wave
fU_current  = fU_in;

% Forward Propagate the Incident Field to Center and NA Crop
if r.useFieldLoss
    fU_inC      = fU_in .* p.propKerH .* p.NA_crop;
end

% Initial Definitions
phsConst    = 2*pi*r.psz/p.lambda;

% Multi-slice Forward Propagation
for idx = 1:size(obj,3)
    arr.eVol(:,:,idx)   = ifft2(fU_current .* p.propKer1);
    fU_current          = fft2(arr.eVol(:,:,idx) .* exp(1i*phsConst*obj(:,:,idx)));
end

% Backward Propagation to Volume Center and NA Crop
efield      = ifft2(fU_current .* conj(p.propKerH) .* p.NA_crop);

% Propagate Field to Defocus Layers and Shift Center Z
propDFZ     = p.propKerDF .* p.propKerZ;
arr.eDF     = ifft2(fft2(efield) .* propDFZ);
if r.useFieldLoss; arr.eInDF = ifft2(fU_inC .* propDFZ); end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                msBwd.m                                  %
% _______________________________________________________________________ %
%              Inverse model helper function for mainMSBP.m               %
%                                                                         %
% This function is a reformulation of BPM_update from previous versions   %
% of the code. It has been reworked to increase execution speed and       %
% improve memory utilization. In addition, the gradient implementation    %
% has been updated to reflect the difference in the real- vs. complex-    %
% valued solutions.                                                       %
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
function [obj, cost] = msBwd(p, r, obj, acq, arr)
%% msBwd: Run inverse gradient update based on multislice beam propagation 
%%        forward model 
% Arguments:
%   p         - System Parameters Struct
%   r         - Reconstruction Parameters Struct
%   obj       - 3D RI volume
%   acq       - Measurement data
%   arr       - Struct holding pre-allocated field volumes
% Returns:
%   obj       - Updated 3D RI volume
%   cost      - Iteration cost update
arguments
    p       struct
    r       struct
    obj     (:,:,:) {mustBeUnderlyingType(obj,'single')}
    acq     (:,:,:) {mustBeUnderlyingType(acq,'single')}
    arr     struct
end

% Calculate Residual Term
if r.useFieldLoss
    backProp    = arr.eDF - acq.*arr.eInDF;
    res         = mean(acq.*arr.eInDF-arr.eDF,3);
else
    backProp    = arr.eDF - abs(acq).*exp(1i*angle(arr.eDF));
    res         = mean(abs(acq)-abs(arr.eDF),3);
end

% Compute Cost
cost        = norm(res(r.pdar+1:end-r.pdar,r.pdar+1:end-r.pdar),'fro')^2;

% Apply Refocusing Kernel and Shift to Center Z
backProp    = ifft2(fft2(backProp) .* conj(p.propKerDF.*p.propKerZ));

% Average Refocused Residual Terms
backProp    = mean(backProp,3);

% Forward Propagate BackProp Term to Final Layer
backProp    = ifft2(fft2(backProp) .* p.propKerH .* p.NA_crop);

% Initial Definitions
phsConst     = 2*pi*r.psz/p.lambda;
conjPropKer1 = conj(p.propKer1);

% Loop
for idx = size(obj,3):-1:1
    % Read object slice
    objSlice    = obj(:,:,idx);

    % Compute gradient
    backProp    = exp(-1i * phsConst * conj(objSlice)) .* backProp;
    grad        = -1i * phsConst * conj(arr.eVol(:,:,idx)) .* backProp;

    % Update object layer
    if r.useComplex
        objSlice    = objSlice - r.stepSize * grad;

        % Apply Imaginary Positivity Constraints (No Gain)
        objSlice    = complex( real(objSlice), max(imag(objSlice),0) );
    else
        objSlice    = objSlice - r.stepSize * real(grad);
    end
    
    obj(:,:,idx)    = objSlice;

    % Propagate residual
    if idx > 1
        backProp = ifft2(conjPropKer1 .* fft2(backProp));
    end
end
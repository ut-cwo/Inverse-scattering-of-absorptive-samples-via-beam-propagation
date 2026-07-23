function [obj, cost] = msBwd(p, r, obj, acq, arr)
%% msBwd: Run backwards gradient update based on multislice model 
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
cost        = norm(res(:))^2;

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

        % Apply Positivity Constraints
        if r.positivityCon
            objSlice    = complex( max(real(objSlice),0),...
                                   max(imag(objSlice),0));
        else
            objSlice    = complex( real(objSlice), ...
                                   max(imag(objSlice),0));
        end
    else
        objSlice    = objSlice - r.stepSize * real(grad);

        % Apply Positivity Constraint
        if r.positivityCon
            objSlice    = max(real(objSlice),0);
        end
    end
    obj(:,:,idx)    = objSlice;

    % Propagate residual
    if idx > 1
        backProp = ifft2(conjPropKer1 .* fft2(backProp));
    end
end
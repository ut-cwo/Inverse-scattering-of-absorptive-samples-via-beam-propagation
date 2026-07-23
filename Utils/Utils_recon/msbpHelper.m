function [reconObj, p, r] = msbpHelper(p, r)
%% msbpHelper: Helper function to package the core msbp loop
% Arguments:
%   p         - System Parameters Struct
%   r         - Reconstruction Parameters Struct
% Returns:
%   reconObj  - Updated 3D RI volume
%   p         - System Parameters Struct
%   r         - Reconstruction Parameters Struct
%% Define physical volume parameters
p.N         = single(p.patchFOV+2*r.pdar);      % Number of pixels in padded FOV
x           = gpuArray(p.ps*(-p.N/2:p.N/2-1));  % Lateral pixel size of object space
[p.xx,p.yy] = meshgrid(x,x);                    % 2D padded grid in x/y
clear x;

p.dfx       = 1/(p.N*p.ps);                     % Lateral pixel size of object frequency space
fx          = gpuArray(p.dfx*(-p.N/2:p.N/2-1)); % 1D padded axis in fx
[p.fxx,p.fyy]   = meshgrid(fx,fx);              % 2D padded grid in fx/fy
clear fx;

%% Propagation Kernel Definition
%   Define frequency support of free space in sample medium and remove
%   evanescent propagation frequencies.
fz              = (p.n_m/p.lambda)^2-(p.fxx.^2+p.fyy.^2);
p.propPhs_m     = 1i*2*pi*sqrt(max(fz,0));
p.propPhs_m     = fftshift(p.propPhs_m);
clear fz;

%% Numerical Aperture Frequency Support
p.NA_crop       = (p.fxx.^2 + p.fyy.^2 < (p.NA/p.lambda)^2);
p.NA_crop       = fftshift(p.NA_crop);

%% Define Propagation Kernels
p.propKer1    = exp(p.propPhs_m .* r.psz);              % Single Layer
p.propKerH    = exp(p.propPhs_m .* r.O/2  .* r.psz);    % Half Volume
p.propKerZ    = exp(p.propPhs_m .* r.zPlane  .* r.psz); % Z Plane Recenter

% Refocus Planes' Kernels
if numel(r.rfDists)==0; r.rfDists = 0; end
p.propKerDF = reshape(single(r.rfDists), 1, 1, []);
p.propKerDF = exp(p.propPhs_m .* p.propKerDF);

%% Initialize Reconstruction Volume
if isempty(r.initGuess)
    reconObj        = zeros([p.N,p.N,r.O],'single','gpuArray');
    if r.useComplex; reconObj = complex(reconObj); end
else
    reconObj    = r.initGuess;
    reconObj    = reconObj(r.pdar+1:end-r.pdar,r.pdar+1:end-r.pdar,:);
    reconObj    = single(gpuArray(reconObj));
    r           = rmfield(r, 'initGuess');
end

proxObj_prev    = reconObj;             % Previous object update for momentum
dualRe          = struct([]);           % Initialize struct to hold real warm start values
dualIm          = struct([]);           % Initialize struct to hold imaginary warm start values
betaDamp        = 0;                    % Initialize beta damping term
t_k             = 1;                    % Initialize tk for nesterov

r.cost          = zeros(r.maxIter,1);   % Monitor cost during convergence

%% Initializing Figure Windows to Observe Iterative Process
close all;

% View reconstruction progress
if r.viewUpdates
    figure('Name','Reconstruction result');
    MSBP_progview(real(reconObj), 1, r.plotRange, r.cost, 0)
    MSBP_progview(imag(reconObj), 2, r.plotRange, r.cost, 0)
    drawnow();
end

%% Define Defocus Planes
% The measured field data is demodulated to DC during processing. To
% correclty refocus, we need to remodulate the signal before propagation.

% Initialize
measAcqs_pad        = padarray(p.measAcqs,[r.pdar,r.pdar,0],1);
propDataset         = complex(zeros(size(measAcqs_pad,1), ...
                                    size(measAcqs_pad,2), ...
                                    size(measAcqs_pad,3), ...
                                    length(r.rfDists), ...
                                         'single','gpuArray'));

% Digitally defocus field measurements
for iii = 1:size(p.measAcqs,3)
    % Define incident field
    fx_in_interp    = fix(p.fx_in(iii)/p.dfx)*p.dfx;
    fy_in_interp    = fix(p.fy_in(iii)/p.dfx)*p.dfx;
    U_in            = exp(1i * 2 * pi * (   fx_in_interp * p.xx ...
                                          + fy_in_interp * p.yy));

    % Modulate incident and measured fields then propagate
    for df = 1:numel(r.rfDists)
        temp1   = ifft2( fft2(measAcqs_pad(:,:,iii).*U_in) ...
                                    .* p.propKerDF(:,:,df));
        temp2   = ifft2( fft2(U_in) .* p.propKerDF(:,:,df));

        % Demodulate data back to DC
        propDataset(:,:,iii,df) = temp1./temp2;
    end
end

% Permute the data as follows: [Y, X, Defocus, Angles]
propDataset = permute(propDataset,[1,2,4,3]);

% For amplitude loss, cast measurements to real
if ~r.useFieldLoss
    propDataset = real(abs(propDataset));
end

% Clear old data
p       = rmfield(p, 'measAcqs');
clear   measAcqs_pad temp1 temp2;

%% MSBP Loop
% Start timer for full reconstruction
r.tLoopCheck = zeros(r.maxIter,1);

% Initialize array to store field volume
arr.eVol    = complex(zeros(p.N,p.N,r.O,'single','gpuArray'));

tTot = tic;
for iter = 1:r.maxIter
    % Start timer for loop
    tLoop = tic;

    if r.randAngles
        % Randomly scramble angles and choose without replacement
        seq = randperm(length(p.fx_in));
    else
        % Span angles sequentially
        seq = 1:(length(p.fx_in));
    end

    for angleScan = 1:length(seq)
        % Skip angles from Omit List
        if ismember(seq(angleScan), r.OmitList)
            continue
        end

        % Create incident planewave
        fx_in_interp    = fix(p.fx_in(seq(angleScan))/p.dfx)*p.dfx;
        fy_in_interp    = fix(p.fy_in(seq(angleScan))/p.dfx)*p.dfx;
        U_in            = exp(1i * 2 * pi * (   fx_in_interp * p.xx ...
                                              + fy_in_interp * p.yy));
        fU_in           = fft2(U_in);

        % Forward simulation of incident wave through current object update
        arr = msFwd(p, r, reconObj, fU_in, arr);
        
        % Compute gradient and update object
        [reconObj, funcVal]  = msBwd(p, r, reconObj, ...
                                     propDataset(:,:,:,seq(angleScan)), ...
                                     arr);

        % Compute accumulated error at current iteration
        r.cost(iter)    = r.cost(iter) + funcVal;

        % Display illuminations completed
        if mod(angleScan,50)==0
            fprintf('Iteration: %1.0d, %1.0d/%1.0d Angles.\n', ...
                                        iter,angleScan,length(p.fx_in))
        end

        % Compute total number of angles scanned
        totAngles  = (iter-1)*length(p.fx_in)+angleScan;

        % Apply acceleration
        if mod(totAngles,r.regFreq)==0
            % Remove first layer object mean
            if r.objMeanSub
                valBot      = mean(real(reconObj(:,:,1)),'all');
                reconObj    = reconObj - valBot;
            end

            % Start timer for regularization
            regT = tic;

            % Proximal TV regularization
            [proxObj, dualRe]   = prox_tv_3d(   real(reconObj), ...
                                                r.regParamRe, ...
                                                r.tvIters, r.warmStart, ...
                                                r.tvItersWarm, dualRe);
            if r.useComplex
                [tmp, dualIm]   = prox_tv_3d(   imag(reconObj), ...
                                                r.regParamIm, ...
                                                r.tvIters, r.warmStart, ...
                                                r.tvItersWarm, dualIm);

                % Recombine real and imaginary regularized components
                proxObj         = complex(proxObj,tmp);
                clear tmp;
            end

            % Print regularization time
            regCheck = toc(regT);
            fprintf('Regularization Completed in %1.3fs\n',regCheck)

            if r.nesterov
                % Nesterov accelerated gradients
                t_k1            = 0.5 * (1 + sqrt(1 + 4 * t_k^2));
                betaNest        = (t_k - 1)/t_k1;
                reconObj        = proxObj + betaNest*(proxObj - proxObj_prev);
                t_k             = t_k1;
                proxObj_prev    = proxObj;
            else
                % Apply momentum with fixed step size
                %   Apply a damping term to lower momentum if cost increases
                betaAcc         = r.beta*min(betaDamp/r.accMom,1);
                reconObj        = proxObj + betaAcc*(proxObj - proxObj_prev);
                proxObj_prev    = proxObj;
                betaDamp        = betaDamp + 1;
                if iter>1
                    if r.cost(iter)>r.cost(iter-1)
                        % Reset momentum
                        betaDamp    = 0;

                        % Extend momentum damping
                        r.accMom    = r.accMom + ceil(r.accMom/10);
                    end
                end
            end
        end
    end

    % Print loop time and error
    r.tLoopCheck(iter) = toc(tLoop);
    fprintf("\nIteration: %d | Completed in: %1.2f (s)\n",iter, r.tLoopCheck(iter))
    fprintf("\tError: %1.2f\n", r.cost(iter));

    costNormFac = p.N^2 * numel(seq);
    fprintf("\tNormalized Pixel Error: %1.4f\n\n", r.cost(iter)/costNormFac);

    if r.viewUpdates
        % Update Display
        MSBP_progview(real(reconObj), 1, r.plotRange, r.cost, iter)
        MSBP_progview(imag(reconObj), 2, r.plotRange, r.cost, iter)
        drawnow();
    end

    % Apply early stopping criterion
    if r.earlyStop && iter>5
        % Check the absolute cost difference of the past 5 iterations
        costDiffSum     = abs(sum(diff(r.cost(iter-4:iter)/costNormFac)));
        fprintf('Cost difference variation between iterations: %1.7f\n\n', costDiffSum)
        if costDiffSum < r.lossThresh
            % Early stopping criteria
            fprintf('Terminated after reaching early stop condition...\n')
            break
        end
    end
end

% Print total time and error
r.tTotCheck = toc(tTot);
fprintf("\nFinal Iteration | Total Time: %1.2f (s)\n", r.tTotCheck)
fprintf("\tAverage Loop Time: %1.2f (s)\n", sum(r.tLoopCheck)/iter)
fprintf("\tFinal Error: %1.2f\n", r.cost(iter));
fprintf("\tFinal Normalized Pixel Error: %1.7f\n\n", r.cost(iter)/costNormFac);

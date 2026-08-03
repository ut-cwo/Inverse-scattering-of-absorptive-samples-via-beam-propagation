%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                             sysParFile.m                                %
% _______________________________________________________________________ %
%                Parameter definition file for mainMSBP.m                 %
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
function r = sysParFile(p)
%% sysParFile:  Function call to define reconstruction parameters
% Arguments:
%   p           - Static dataset parameters
% Results:
%   r           - Reconstruction parameters
% _______________________________________________________________________ %
%                                                                         %
%                    [General Reconstruction Variables]                   %
% _______________________________________________________________________ %
r.rfDists   =[-10 0 10];    % Refocus distances [um]
% Define the number of defocus planes and the amount of defocus
% e.g.  [-10 0 10] would be three defocus planes at -10um, 0um, and 10um 
% from the focal plane.
% Larger distances are observed to be needed for thin, highly absorptive
% samples.
%   Common values:    [-10 0 10], [-5 0 5]

r.maxIter   = 200;          % Maximum number of iterations
r.stepSize  = 1e-5;         % Gradient step size
% Larger values tend to cause the solution to diverge for complex-valued 
% reconstructions.
%   Common values:    1e-5

r.OmitList  = [];           % List of measurements to remove
                            %   (Dust, artefacts, etc.)

% _______________________________________________________________________ %
%                                                                         %
%                   [Reconstruction Volume Parameters]                    %
% _______________________________________________________________________ %
r.O         = 200;      % Number of layers for reconstruction
r.psz       = 4*p.ps;   % Size of axial layers
% Number of layers and psz define the axial volume size (psz*O). Lower O
% values means the reconstruction will be faster but for the same volume
% thickness, psz will need to increase. This leads to less axial resolution 
% and can cause reconstruction artefacts as O gets too small.
%   Common values:    psz = 3*ps, 5*ps

r.zPlane    = 25;       % Define center plane of reconstruction volume
% If the measured images are out of focus from the center of the object, 
% the volume reconstruction will not be centered. Instead of increasing O
% and hence, the compute time, you can shift the volume towards the center.
%   Positive values shifts the object up by 'zPlane' layers
%   Negative values shifts the object down by by 'zPlane' layers

r.pdar      = 0;        % Padding size to avoid edge artifacts
%   Choice of padding is important:
%       Large padding -> slow compute or OOM
%       Small padding -> boundary artefacts

r.initGuess = [];       % Initial volume guess
% Define a non-zero intial guess for the reconstructed object. This could 
% be a previous reconstruction or coarse recosntruction. Confirm the grid
% and pixel sizes of the loaded object match the reconstruction parameters.

% _______________________________________________________________________ %
%                                                                         %
%                 [Real- vs. Complex-valued Reconstruction]               %
% _______________________________________________________________________ %
% Define wether to reconstruct a real-valued refractive index volume or a
% complex-valued refractive index volume.
r.useComplex    = true;   
%   True:     Reconstruct complex-valued gradients.
%   False:    Reconstruct real-valued gradients.

% _______________________________________________________________________ %
%                                                                         %
%                    [Amplitude- vs. Field-based loss]                    %
% _______________________________________________________________________ %
% Use an amplitude- or field-based loss for the gradient optimization.
% Field-based loss measures the difference between the measured and 
% simulated fields. Hence, it requires a field-based measurement system. It 
% has been observed, however, that for thick heterogenous samples, field-
% based loss performs worse than ampltiude-based loss due to phase wrapping
% artefacts. This code is designed to be used with field-based measurements
% but the loss can be performed with either amplitude or field loss. Due to
% the presence of field, the amplitude measurements can be digitally
% refocused to encode phase information into the amplitude gradients to
% improve the reconstruction without using a field-based loss.
r.useFieldLoss  = false;
%   True:     Use field-based loss.
%   False:    Use amplitude-based loss.

% _______________________________________________________________________ %
%                                                                         %
%                [Regularization and Constraint Parameters]               %
% _______________________________________________________________________ %
r.regParamRe    = 10e-5;     % Regularization strength RI
r.regParamIm    = 10e-5;     % Regularization strength Absorptivity
% Regularization strength for real and imagnary (RI and Absorption)
% components. Can be tuned separately.
% Regularization strengths are tied to the step size and volume size as 
% they are not normalized. Tuning the step size up or down or changing the 
% patch size may require tuning the regularization stengths.
%   Common values:    10e-5, 20e-5.

r.regFreq       = 400;       % How often to regularize (in angles)
% The regularization frequency term defines how often to regularize and
% apply boundary constraints. Commonly, this is to set to be equal to the
% number of angles in the dataset so regularization is applied each 
% iteration.

r.tvIters       = 100;       % Number of standard proximal TV iterations
%   Common value:     100

r.warmStart     = true;
% Apply a warm start criteria to the proximal TV regularization. This
% stores the previous dual values to increase the speed of proximal
% convergence. Increases memory usage but decreases regularization time.
%   True:     Use warm-start for the proximal TV regularization.
%   False:    Use standard TV.
r.tvItersWarm   = 20;        % Number of TV iterations with warm start
% Define how many iterations to apply for the proximal TV step if using the
% warm-startint criteria. Because the warm start is closer to a good
% proximal value, this can be much smaller than the standard TV iterations.
%   Common value:     10

r.positivityCon = false;
% Positivty constraint can be useful for making the volume pop-out from the
% background. Can also help with convergence.
%   True:     Apply positivity to the Refractive Index reconstruction.
r.objMeanSub    = false;  
% Mean value subraction removes the real-valued refractive index background
% by subtracting the mean value of the first layer from the whole object.
%   True:     Remove mean from object volume during reconstruction.

% _______________________________________________________________________ %
%                                                                         %
%                     [Accelerated Gradient Parameters]                   %
% _______________________________________________________________________ %
r.nesterov      = false;
% True:     Use nesterov acceleration
% False:    Use standard momentum
r.beta          = .999;      % Standard momentum strength
%   Common values: .9,.99,.999
r.accMom        = 5;        % Momentum acceleration dampling
% Ramp up momentum gradually across 'accMom' steps. If cost increases, beta
% is reset and accMom is incremented.
%   Common value: 5

% _______________________________________________________________________ %
%                                                                         %
%                           [Stopping Criterion]                          %
% _______________________________________________________________________ %
% A check to see if the absolute cost difference between the last 5 
% iterations is below some threshold. When true, terminate the
% optimization loop early.
r.earlyStop     = true;
% True: Enable early stopping criterion for code under loss threshold
r.lossThresh    = 1e-5;     % Absolute difference threshold
%   Common values: 1e-5, 1e-6

% _______________________________________________________________________ %
%                                                                         %
%                     [Randomize sequence of gradients]                   %
% _______________________________________________________________________ %
% Random angle sequencing during optimization may increase convergence
% speeds.
r.randAngles    = false;  
%   True:     Iterate over angles randomly during update.
%   False:    Iterate over angles sequentially (1-N) during updates.

% _______________________________________________________________________ %
%                                                                         %
%                             [Visualization]                             %
% _______________________________________________________________________ %
% Allow the user to view the progress of the reconstruction each iteration.
r.viewUpdates   = true;
%   True:     View reconstruction after each iteration
r.plotRangeRe   = [-0.005,0.04];    % Colorbounds for real-valued RI
r.plotRangeIm   = [0,0.02];         % Colorbounds for complex-valued RI

end

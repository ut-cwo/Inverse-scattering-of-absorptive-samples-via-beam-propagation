%% Main reconstruction code
% Main script for running MSBP on measured field data for complex RI 
% reconstructions. Field measurements can be used with either an amplitude- 
% or field-based loss function. For this work, we use digitally
% defocused amplitude measurements from field data for the optimization to 
% move away from field-based loss due to artefacts we observe when using 
% field data in thick samples [1,3].
%
% This code is a reworked version of [2-3] implementing the necessary 
% changes for absorption based reconstruction. Notable changes to the 
% orignal code structure were also made to improve reconstructions speeds. 
% This code utilizes regularization using a 3D total variation (TV) 
% proximal operator. We gratefully acknowledge the work done by Beck and 
% Teboulle [4] for the derivation of proximal tv in 3D. The proximal TV 
% code [4] used in previous works [2,3] was modified for compute speed 
% purposes by adding a warm start capability and other small adjustments.
%
% If using this code, please cite [1].
%
% Author: Peter Wagenaar, Jeongsoo Kim, Shwetadwip Chowdhury; July 22, 2026
%
% References for MSBP:
% [1]   Wagenaar, Peter, et al. "Inverse-scattering of absorptive samples 
%       via beam propagation." bioRxiv (2026).
%
% [2]   S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and L. 
%       Waller, "High-resolution 3D refractive index microscopy of 
%       multiple-scattering samples from intensity images," 
%       Optica 6, 1211-1219 (2019)
%
% [3]   Kim, Jeongsoo, et al. "Inverse-scattering in biological samples via 
%       beam-propagation." bioRxiv (2025).
%
% [4]   A. Beck and M. Teboulle. Fast gradient-based algorithms for 
%       constrained  total variation image denoising and deblurring 
%       problems. Image Processing, IEEE Transactions on, 
%       18(11):2419--2434, 2009
%
% Acknowledgements:
% Thank you to Michael Chen and David Ren for preliminary versions of this
% code.

%% Instructions
% In windows, MATLAB can be run block by block using ctrl + enter or the
% run and advance button in the 'SECTION' tab under 'EDITOR'
%
% This code has been developed exclusivley with gpu computing in mind. 
% Hence, a GPU and the Parallel Computing Toolbox are required to run this
% code. You can download the toolbox by going to the 'HOME' tab and
% selecting 'Add-Ons'. Once there, search for 'Parallel Computing Toolbox'
% and install.
%
% To prepare the code, first follow the steps below:
%   1.  Download the data from the Texas Data Repository and move the code
%       to a new folder. MATLAB will run within this folder and all
%       references will reference this path or a fullfile path.
%       Code structure should be as below:
%           ...\Utils\
%               ...\Utils_IO
%               ...\Utils_vis
%               ...\Utils_reg
%               ...\Utils_recon
%               ...\Utils_ffts
%           ...\mainMSBP.m
%           ...\msbpHelper.m
% 
%   2.  Next, set the 'p.dataPath' variable in [Block 1] to the location of 
%       your data.
%           e.g.    p.dataPath = 'D:\User\Data\Zebrafish\' or 'Zebrafish' 
%                   if the folder is located on the MATLAB path
%       This does not need to be in the same location as the code as long
%       as you reference the full path.
%
%   3.  [Block 1] also has the reference to the sample name. Set p.sampName
%       to the name of the path of the field measurements.
%           e.g.    p.sampPath = 'sample_fields_FOV1.mat'
%       
%   4.  Next, set the reconstruction parameters. This can be done by 
%       loading the accompanying 'parFile.mat' or by manually setting the
%       variables. To load parameters from the 'parFile', just set the path
%       of the parFile variable in [Block 2] to the associated
%       'parFile.mat'
%           e.g.    parFile = 'D:\User\Data\Zebrafish\parFile_ZF.mat' or 
%                   'parFile_ZF.mat' if the folder is located on the MATLAB 
%                   path
%       This does not need to be in the same location as the code as long
%       as you reference the full path.
%
%       To manually set the variables, change the values of the
%       corresponding value in r.
%           e.g.    r.O = 200; r.useComplex = false; etc.
%
%   5.  Next, define the patch size and position for reconstrution. In
%       [Block 3], set 'p.patchFOV', 'p.xCent', and 'p.yCent' such that the
%       region of interest is highlighted by the red box in the image.
%       Values can be tuned in this box and checked by continously
%       rerunning the same block.
%
%   6.  Finally, run the msbpHelper function in [Block 4]. This will read
%       the data ('p') and reconstruction ('r') parameters and run the
%       multislice code. The function outputs the updated structs as well
%       as the reconstructed object.
%
%   7. To save the data, run [Block 5]. This will create a new folder
%       storing the reconstructed object, parameters, metadata, and code
%       files following [Block 5].
%
%   Notes:  A visualization tool 'sliderDisplayImVC2' is included to easily
%           view the 3D measurements and reconstructed object.
%               e.g. sliderDisplayImVC2(data); colormap gray; clim([0,1.5])
%                       data: 3D volume
%
%   To run this code with external datasets, a .mat file should be created
%   in the following format:
%       Efield_amplitude    -   3D numeric array holding amplitude for each 
%                               measurement. [Y X Angles] (float)
%       Efield_phase        -   3D numeric array holding phase for each
%                               measurement. [Y X Angles] (float)
%       fx_illum_ref        -   2D vector holding fx components of 
%                               reference angles. [1 Angles] (float)
%       fy_illum_ref        -   2D vector holding fy components of 
%                               reference angles. [1 Angles] (float)
%       ps                  -   Pixel size in [um] (float)
%       lambda              -   Wavelength in [um] (float)
%       NA                  -   Numerical aperture (float)
%       n_m                 -   Refractive index of sample media (float)
%       n_imm               -   Refractive index of immersion (float)
%
%   Field measurements from the accompanying datasets contain background
%   subtracted field measurements with accompanying spatial frequency
%   components. Imported data shuould follow this convention.

%% Clear MATLAB Workspace
clear
clc
close all

%% Add Relevant Paths
addpath(genpath('Utils'))

%% [Block 1] Define Data Paths [Requires User Definitions]
p.dataPath  = '';
p.sampName  = '';

%% Load in Data
% Load in field measurements
data        = loadData(p.dataPath, p.sampName);

% When loading external measurements, you may find it necessary to apply a 
% data correction term to flip the phase of the field measurements. This 
% comes from angle calibration inconsistencies. If the reconstructed RI 
% looks inverted, try setting flip phase to -1.
p.flipPhase = 1;

% Combine amplitude and phase to get complex field measurements.
totFOV_acqs = data.Efield_amplitude.*exp(1i.*p.flipPhase*data.Efield_phase);

% Load K vectors
%   Similar to above, if reconstruct looks flipped or is poorlt
%   reconstructed, the calibrated kx, ky values may need to be flipped.
p.fx_in     = data.fx_illum_ref;
p.fy_in     = data.fy_illum_ref;

% Load physical parameters
p.ps        = data.ps;
p.lambda    = data.lambda;
p.NA        = data.NA;
p.n_m       = data.n_m;
p.n_imm     = data.n_imm;

% Display System Variables
fprintf( ...
    ['------------------\n' ...
     'ps = %2.3f\n' ...
     'lambda = %2.3f\n' ...
     'NA = %2.2f\n' ...
     'n_m = %2.3f\n' ...
     'n_imm = %2.3f\n' ...
     '------------------\n'], p.ps, p.lambda, p.NA, p.n_m, p.n_imm)

%% [Block 2] Define Object Reconstruction Parameters [Requires User Definitions]
% Read in Parameters from File
parFile       = '';   
% Leave parFile empty for manually inputed recosntruction variables.

if ~isempty(parFile)
    load(parFile)
else
    %% Reconstruction Variables
    r.rfDists   =[-5 0 5];% Refocus distances [um]
                            % Define the number of defocus planes and the 
                            % amount of defocus
                            % e.g.  [-10 0 10] would be three defocus 
                            %       planes at -10um, 0um, and 10um from the 
                            %       focal plane.

    r.O         = 200;      % Number of layers for reconstruction
    r.psz       = 6*p.ps;   % Size of axial layers
  
    r.zPlane    = zPlane;   % Shift center plane of reconstruction volume
                            %       Positive values shift object up
                            %       Shifts by X number of planes
    
    r.maxIter   = 50;       % Maximum number of iterations
    r.stepSize  = 1e-5;     % Gradient step size (1e-5 common)

    r.OmitList  = [];       % List of measurements to remove (Dust, artefacts, etc.)


    %% Padding
    r.pdar      = 100;      % Padding size to avoid edge artifacts
                            %   Choice of padding is important:
                            %       large padding -> slow compute or OOM
                            %       small padding -> boundary artefacts

    %% Init
    r.initGuess = [];       % Define a non-zero intial guess for the 
                            % reconstructed object.

    %% Real- vs. Complex-valued Reconstruction
    r.useComplex    = true;
    % True:     Reconstruct complex-valued gradients.
    % False:    Reconstruct real-valued gradients.
    
    %% Amplitude- vs. Field-based loss
    r.useFieldLoss  = false;  
    % True:     Use field-based loss.
    % False:    Use amplitude-based loss.
    
    %% Boundary Conditions
    r.positivityCon = true;
    % True:     Apply positivity to the Refractive Index reconstrucion.
    r.objMeanSub    = false;  
    % True:     Remove mean of first layer from object update during reconstruction.
    
    %% Regularization Parameters
    r.regParamRe    = 10e-5;     % Regularization strength RI
    r.regParamIm    = 10e-5;     % Regularization strength Absorptivity
    % Regularization strength for real and imagnary (RI and Absorption)
    % components. Can be tuned separately. 
    r.regFreq       = 100;       % How often to regularize (in angles)
    r.tvIters       = 100;       % Number of TV iterations w/o warm start
    r.tvItersWarm   = 20;        % Number of TV iterations after warm start
    r.warmStart     = true;
    % True:     Uses an initial guess for the proximal TV regularization.
    %               Increases memory usage but speeds up reconstruction.
    
    %% Accelerated Gradient Parameters
    r.nesterov      = false;
    % True:     Use nesterov acceleration
    % False:    Use standard momentum
    r.beta          = .999;      % Standard momentum strength
    r.accMom        = 10;       % Ramp up momentum gradually across accMom
                                %   If cost starts to increase, momentum
                                %   is reset and accMom grows.
    
    %% Stopping Criterion
    r.earlyStop     = true;
    % True: Enable early stopping criterion for code under loss threshold
    r.lossThresh    = 1e-5; % Percent difference threshold for previous costs
    
    %% Randomize Sequence of Gradients
    r.randAngles    = false;  
    % True:     Use random sequence of angles during update.
    %               Random anglular sequencing may speed up convergence.
    % False:    Use sequential angulese during updates.

    %% Visualization
    r.plotRange     = [-0.01,0.05]; % Colorbounds for progress update
    r.viewUpdates   = true;
    % True:     View reconstruction after each iteration

    %% Save Paramater Struct as .mat
    if ~exist([p.dataPath 'ParFiles/'],'dir')
        mkdir([p.dataPath 'ParFiles/'])
    end
    save([p.dataPath 'ParFiles/parFile.mat'],'r');
end

%% [Block 3] Define Patch FOV [Requires User Definitions]
% Define the patch you want to reconstruct
p.patchFOV      = 0;     % Define patch size
p.xCent         = 0;      % Define x center
p.yCent         = 0;      % Define y center

% Assert sizes are even
assert(mod(p.patchFOV,2)==0,'Transverse Patch Size should be even...')
assert(mod(r.O,2)==0,'Layers should be even...')

% Define rows and coloumns of patch
rows            = p.yCent+(-p.patchFOV/2+1:p.patchFOV/2);
cols            = p.xCent+(-p.patchFOV/2+1:p.patchFOV/2);

% Plot
figure;
imagesc(abs(totFOV_acqs(:,:,1))); clim([0,2]); 
axis equal; axis tight; 
colormap gray; 
hold on;
rectangle('Position',[cols(1),rows(1),p.patchFOV,p.patchFOV], ...
          'LineWidth',3, 'EdgeColor','r'); 
hold off;

p.measAcqs  = gpuArray(totFOV_acqs(rows,cols,:));

%% Optimal patch FOV with padding
% To maximize speed of compute, choosing the correct size of computed
% fields can make a big difference. As multi-slice beam propagation heavily 
% relies on FFTs, running them inefficiently can greatly affect the compute 
% time.
%
% Optimal cuFFT size requirements:
%   Size in form of 2^a*3^b*5^c*7^d
%       In general the smaller the prime factor, the better the performance.
%           Powers of two will be the fastest.

fastPrimes  = [2 3 5 7];
fovFactors  = factor(p.patchFOV+2*r.pdar);

fprintf('Patch FOV with padding: %d\n', size(p.measAcqs,1)+2*r.pdar)
fprintf('Factorization: \t\t\t'); 
fprintf('%d ', fovFactors);
fprintf('\n')

if ~all(ismember(fovFactors, 2),'all')
    warning("Size is not a power of 2. " + ...
            "Will not achieve fastest compute time...")
end
if ~all(ismember(fovFactors, fastPrimes),'all')
    warning("Size is not a factorization of only 2, 3, 5, and 7. " + ...
            "Compute times may be much slower...")
end

%% Note: Be mindful of memory utilization
% Proximal TV is memory hungry and can transiently exceed the total VRAM on
% GPU. This can throw an OOM error or silently cause the reconstruction to 
% slow down by writing to and from shared memory.
% You can check if this is happening in the task manager from the dedicated
% and shared memory utilization.
% If you notice the GPU writing to shared memory and want to increase 
% reconstruction speed, you can patch sample into smaller FOVs and run 
% sequentially. This can be faster than running too large of an FOV if
% implemented efficiently.

%% Clear After Patch Found
clear data totFOV_acqs rows cols;

%% Define Path to Save and Variables to Save
reconsDir = p.dataPath;
r.varsToSave = {'p', 'r'};

%% [Block 4] Run Loop
% Arguments:
%   p:  Struct with all the physical system parameters, these remain
%       constant between runs.
%   r:  Struct containing reconstruction variables. These will change
%       under different reocnstruction params.
% Returns:
%   reconObj:   Reconstructed complex/real-valued RI map.
%   p:          Struct with all the physical system parameters, these 
%               remain constant between runs.
%   r:          Updated struct containing reconstruction variables. These 
%               will change under different reocnstruction params.

[reconObj, p, r] = msbpHelper(p, r);

%% [Block 5] Save
disp(['Saving data to : ' reconsDir]);

saveRecon(reconsDir, reconObj, r.varsToSave, ...
    'Format','mat', 'IDMode','number', ...
    'CodeFiles', {'msFwd','msBwd','msbpHelper'}, ...
    'Metadata', struct('Notes',''));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               mainMSBP.m                                %
% _______________________________________________________________________ %
%                 Main script used to launch msbpHelper.m                 %
%                                                                         %
% Main script for running the multi-slice beam propgation inverse solver  %
% using measured field data for complex-valued RI reconstructions. In     %
% this work, field measurements are used to create computationally        %
% defocused amplitude measurements. This is chosen over using a field-    %
% based loss to remove artefacts observed from field-based loss in thick  %
% samples [1,2].                                                          %
%                                                                         %
% This code is a reworked version of [2-3] implementing the necessary     %
% changes for absorption based reconstruction. In addition, notable       %
% changes to the orignal code were made to improve reconstruction speeds. % 
%                                                                         %
% Finally, this code incorporates regularization using a 3D total         %
% variation (TV) proximal operator. We gratefully acknowledge the work    %
% done by Beck and Teboulle [4] for the derivation of proximal tv in 3D.  %
% The proximal TV code used in previous works [2,3] was modified for      %
% compute speed purposes by adding a warm start capability and increasing %
% efficiency for large gpuArrays.                                         %
%                                                                         %
% If using this code, please cite [1].                                    %
%                                                                         %
%                                                                         %
% Authors:                                                                %
%   Peter Wagenaar, Jeongsoo Kim, Shwetadwip Chowdhury; July 22, 2026     %
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
%                                                                         %
%                                                                         %
% [4]   A. Beck and M. Teboulle. Fast gradient-based algorithms for       %
%       constrained  total variation image denoising and deblurring       %
%       problems. Image Processing, IEEE Transactions on,                 %
%       18(11):2419--2434, 2009                                           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                         Instructions for the User                       %
% ________________________________________________________________________%
% In windows, MATLAB can be run section by section with "ctrl + enter"    %
% while focused on the corresponding section. These sections are          %
% defined between instances of '%%'. They can also be run using the       %
% "Run and Advance" button in the 'SECTION' tab under 'EDITOR'.           %
%                                                                         %
% ** IMPORTANT**                                                          %
% This code has been developed exclusivley with GPU computing in mind.    %
% Hence, both a GPU and the Parallel Computing Toolbox are required. You  %
% can download the toolbox by going to the 'HOME' tab and selecting       %
% 'Add-Ons'. Once there, search for 'Parallel Computing Toolbox' and      %
% install.                                                                %
%                                                                         %
% Confirm the code has been dowloaded correctly by matching the following %
% folder structure:                                                       %
%           ...\Utils\                                                    %
%                   \Utils_IO                                             %
%                   \Utils_recon                                          %
%                   \Utils_reg                                            %
%                   \Utils_vis                                            %
%           ...\ParameterFiles\                                           %
%           ...\mainMSBP.m                                                %
%                                                                         %
% Next, to prepare to run this code, please follow the steps below:       %
%                                                                         %
%   1.  Download the expiremental data from the Texas Data Repository     %
%       link found in the GitHub repository README. This can be placed    %
%       anywhere on your system as long as the full path is referenced.   %
%       Anything placed directly in the same folder as mainMSBP.m will be %
%       added to the path and can be referenced locally.                  %
%                                                                         %
%                                                                         %
%   2.  Set the 'p.dataPath' variable in [Block 2] to the location of     %
%       your data.                                                        %
%                                                                         %
%       e.g.    p.dataPath  = 'D:\User\Data\ZebrafishEmbryo\'             %
%                                   or                                    %
%                           = 'ZebrafishEmbryo\' (Local)                  %
%                           = '' (Local)                                  %
%                                                                         %
%                                                                         %
%   3.  Set the 'p.sampName' variable in [Block 2] to the name of the     %
%       sample.                                                           %
%                                                                         %     
%       e.g.    p.sampPath = 'sampleField_ZF_FOV1.mat'                    %
%                                                                         %
%                                                                         %
%   4.  Set the reconstruction parameters. This can be done by specifying %
%       a path in 'parFile' from [Block 3] or by leaving 'parFile' empty  %
%       ('') and setting the variables in the accompanying 'sysParFile.m' %
%       file located under '\ParameterFiles\'.                            %
%                                                                         %
%       e.g.    parFile = 'D:\User\Data\Zebrafish\parFile_ZF.mat'         %
%                                   or                                    %
%                       = 'parFile_ZF.mat' (Local)                        %
%                                                                         %
%       To manually set the variables:                                    %
%               parFile = '';                                             %
%                                                                         %
%               \ParameterFiles\sysParFile                                %
%                           r.O = 200;                                    %
%                           r.maxiter = 100;                              %
%                                                                         %
%                                                                         %
%   5.  Define the patch size and position for reconstruction in the FOV. %
%       In [Block 4], set 'r.patchFOV', 'r.xCent', and 'r.yCent' and run  %
%       the section. You will see a red box highlight the choosen patched %
%       region. To change this position, alter the values accordingly and %
%       rerun the section until the intended patch is highlighted.        %
%                                                                         %
%                                                                         %
%   6.  Run the section indicated by [Block 9]. This will input the data  %
%       ('p') and reconstruction ('r') parameters into the msbpHelper     %
%       function and run the multi-slice code. The function outputs the   %
%       updated structs as well as the reconstructed object.              %
%                                                                         %
%                                                                         %
%   7.  To save the data, run the section indicated by [Block 10]. This   %
%       will create a new folder storing the reconstructed object,        %
%       reconstruction and data structs, parameter file, metadata, and    %
%       code files.                                                       %
%                                                                         %
%                                                                         %
%   Notes:  A visualization tool 'sliderDisplayImVC2' is included to      % 
%           easily view the 3D measurements and reconstructed object.     %
%               e.g.    sliderDisplayImVC2(data);                         %
%                           colormap gray; clim([0,1.5]);                 %
%                       data: 2D angular measurements, 3D volume          %
%                                                                         %
%   To run this code with external datasets, a .mat file should be        %
%   created in the following format:                                      %
%       Efield_amplitude    -   3D numeric array holding amplitude for    % 
%                               each measurement. [Y X Angles] (float)    %
%       Efield_phase        -   3D numeric array holding phase for each   %
%                               measurement. [Y X Angles] (float)         %
%       fx_illum_ref        -   2D vector holding fx components of        %
%                               reference angles. [1 Angles] (float)      %
%       fy_illum_ref        -   2D vector holding fy components of        %
%                               reference angles. [1 Angles] (float)      %
%       ps                  -   Pixel size in [um] (float)                %
%       lambda              -   Wavelength in [um] (float)                %
%       NA                  -   Numerical aperture (float)                %
%       n_m                 -   Refractive index of sample media (float)  %
%       n_imm               -   Refractive index of immersion (float)     %
%                                                                         %
%   Field measurements from the accompanying datasets contain background  %
%   subtracted field measurements with accompanying spatial frequency     %
%   components. Imported data shuould follow this convention.             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% _____________________________________________________________________ %%
%%                    [Block 1] Initialize Workspace                     %%
%% _____________________________________________________________________ %%
% Clear MATLAB Workspace
clear
clc
close all

% Add Relevant Paths
addpath(genpath('Utils'))
addpath(genpath('ParameterFiles'))

%% _____________________________________________________________________ %%
%%                 [Block 2] Define Paths and Load Data                  %%
%%                      [Requires User Definitions]                      %%
%% _____________________________________________________________________ %%
p.dataPath  = '';   % [REQUIRES USER INPUT]
p.sampName  = '';   % [REQUIRES USER INPUT]

% Load in field measurements
fprintf("Loading data from file '%s%s'...\n",p.dataPath, p.sampName)
data        = loadData(p.dataPath, p.sampName);
fprintf("Finished loading data.\n")

% When loading external measurements, you may find it necessary to apply a 
% data correction term to flip the phase of the field measurements. This 
% comes from angle calibration inconsistencies. If the reconstructed RI 
% looks inverted, try setting flip phase to -1.
p.flipPhase = 1;

% Combine amplitude and phase to get complex field measurements.
totFOV_acqs = data.Efield_amplitude.*exp(1i.*p.flipPhase*data.Efield_phase);

% Load K vectors
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

%% _____________________________________________________________________ %%
%%           [Block 3] Define Object Reconstruction Parameters           %%
%%                      [Requires User Definitions]                      %%
%% _____________________________________________________________________ %%
% Read in Parameters from File
parFile       = '';

% Leave parFile empty for manually inputed reconstruction variables.
if ~isempty(parFile)
    fprintf('Loading external paramater file: %s...\n', parFile);
    r = load(parFile);
    fprintf('Finished loading reconstruction parameters:\n');
    disp(r);
else
    fprintf('Loading system paramater file...\n');
    r = sysParFile(p);
    fprintf('Finished loading reconstruction parameters:\n');
    disp(r);
end

%% _____________________________________________________________________ %%
%%                      [Block 4] Define Patch FOV                       %%
%%                      [Requires User Definitions]                      %%
%% _____________________________________________________________________ %%
% Define the patch you want to reconstruct
        % [REQUIRES USER INPUT] %
r.patchFOV      = ;         % Define patch size
r.xCent         = ;         % Define x center
r.yCent         = ;         % Define y center
        % [REQUIRES USER INPUT] %

% Assert sizes are even
assert(mod(r.patchFOV,2)==0, 'Transverse Patch Size should be even...')
assert(mod(r.O,2)==0,        'Layers should be even...')

% Define rows and coloumns of patch
rows            = r.yCent+(-r.patchFOV/2+1:r.patchFOV/2);
cols            = r.xCent+(-r.patchFOV/2+1:r.patchFOV/2);

% Plot
figure;
imagesc(abs(totFOV_acqs(:,:,1))); 
clim([0,2]); axis equal; axis tight; colormap gray; hold on;
rectangle( 'Position',[cols(1),rows(1),r.patchFOV,r.patchFOV], ...
           'LineWidth',3, 'EdgeColor','r'); 
hold off;

p.measAcqs  = gpuArray(totFOV_acqs(rows,cols,:));

%% _____________________________________________________________________ %%
%%            [Block 5] Note: Optimal patch FOV with padding             %%
%% _____________________________________________________________________ %%
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
fovFactors  = factor(r.patchFOV+2*r.pdar);

fprintf('Patch FOV with padding: %d\n', size(p.measAcqs,1)+2*r.pdar)
fprintf('Factorization: \t\t\t'); 
fprintf('%d ', fovFactors);
fprintf('\n')

if ~all(ismember(fovFactors, 2),'all')
    warning("Size is not a power of 2.")
end
if ~all(ismember(fovFactors, fastPrimes),'all')
    warning("Size is not a factorization of only 2, 3, 5, and 7.")
end

%% _____________________________________________________________________ %%
%%           [Block 6] Note: Be mindful of memory utilization            %%
%% _____________________________________________________________________ %%
% Proximal TV is memory hungry and can transiently exceed the total VRAM on
% GPU. This can throw an OOM error or silently cause the reconstruction to 
% slow down due to writing between shared and dedicated GPU memory.
% You can check if this is happening in the task manager by monitoring if 
% there is any shared memory utilization.
% If this occurs during reconstruction, it is highly recommended to use a
% smaller patch or fewer layers to improves reconstruction speeds.

%% _____________________________________________________________________ %%
%%           [Block 7] Clean Up Workspace After Patch is Found           %%
%% _____________________________________________________________________ %%
clear data totFOV_acqs rows cols;

%% _____________________________________________________________________ %%
%%              [Block 8] Define Path and Variables to Save              %%
%% _____________________________________________________________________ %%
reconsDir = [p.dataPath ''];

% Make directory
if ~isempty(reconsDir) && ~exist(reconsDir,'dir')
    mkdir(reconsDir)
end

fprintf('Save directory set: %s\n', reconsDir);

r.varsToSave = {'p', 'r'};

%% _____________________________________________________________________ %%
%%                          [Block 9] Run Loop                           %%
%% _____________________________________________________________________ %%
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

%% _____________________________________________________________________ %%
%%                            [Block 10] Save                            %%
%% _____________________________________________________________________ %%
disp(['Saving data to : ' reconsDir]);

saveRecon(reconsDir, reconObj, r.varsToSave, ...
    'Format','mat', 'IDMode','number', ...
    'CodeFiles', {'msFwd','msBwd','msbpHelper'}, ...
    'Metadata', struct('Notes',''));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           MSBP_progview.m                               %
% _______________________________________________________________________ %
%                   Progress viewer for MSBP updates                      %
%                                                                         %
% Displays a triframe cross-sectional view of a 3D volume, and may also   %
% plot a figure associated with a cost curve (optional)                   %
%                                                                         %
%                                    July 6, 2026 by Shwetadwip Chowdhury %
%                                                                         %
% References:                                                             %
% [3]   S. Chowdhury, M. Chen, R. Eckert, D. Ren, F. Wu, N. Repina, and   %
%       L. Waller, "High-resolution 3D refractive index microscopy of     %
%       multiple-scattering samples from intensity images,"               %
%       Optica 6, 1211-1219 (2019)                                        %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/Waller-Lab/multi-slice                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MSBP_progview(obj,fignum,plotRange,pltCurve,upto)
%% MSBP_progview: View iterative progression
% Arguments:   
%   obj         - 3D matrix 
%   fignum      - Number assigned to this figure
%   plotRange   - Color axis with which to display values
%   pltCurve    - Curve that may get plot (optional)
%   upto        - Boundary up to which curve will be plot
figure(fignum);
subplot(2,2,1);
imagesc(real(squeeze(obj(:,:,end/2)))); axis equal; axis tight;
clim(plotRange); colormap gray; colorbar; title('x,y');
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,2,2);
imagesc(real(squeeze(obj(:,end/2,:)))); axis equal; axis tight;
clim(plotRange); colormap gray; colorbar; title('x,z');
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,2,3);
imagesc(real(squeeze(obj(end/2,:,:)))); axis equal; axis tight;
clim(plotRange); colormap gray; colorbar; title('y,z');
set(gca,'xtick',[])
set(gca,'ytick',[])

subplot(2,2,4);
plot(log10(1+pltCurve),':o'); title('cost function'); axis tight;

if nargin == 5
    subplot(2,2,4);
    start = 1;
    plot(log10(1+pltCurve(start:upto)),':o'); title('cost function'); axis tight;
end
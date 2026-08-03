%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grad_3D.m                                 %
% _______________________________________________________________________ %
%                3D gradient implementation using circshift               %
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
function [dx, dy, dz] = grad_3d(x)
%% grad_3d: 3D gradient using circshift
% Arguments:
%   x       - 3D volume
% Returns:
%   dx      - 3D gradient along x
%   dy      - 3D gradient along y
%   dz      - 3D gradient along z
dx = circshift(x,-1,1) - x;  dx(end,:,:) = 0;
dy = circshift(x,-1,2) - x;  dy(:,end,:) = 0;
dz = circshift(x,-1,3) - x;  dz(:,:,end) = 0;
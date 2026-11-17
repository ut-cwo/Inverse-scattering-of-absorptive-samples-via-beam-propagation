%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                div_3D.m                                 %
% _______________________________________________________________________ %
%               3D divergence implementation using circshift              %
%                                                                         %
%                                                                         %
%                                         July 22, 2026 by Peter Wagenaar %
%                                                                         %
% References:                                                             %
% [1]    Peter Wagenaar, Jeongsoo Kim, Mary E. Swartz,                    %
%        Johann K. Eberhart, and Shwetadwip Chowdhury, "Inverse           %
%        scattering of absorptive samples via beam propagation,"          %
%        Optica 13, 1891-1902 (2026)                                      %
%                                                                         %
%       (Github Repository):                                              %
%           https://github.com/ut-cwo/Inverse-scattering-in-biological-   %
%               samples-via-beam-propagation                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = div_3d(px, py, pz)
%% div_3d: 3D divergence using circshift
% Arguments:
%   px      - 3D field along x
%   py      - 3D field along y
%   pz      - 3D field along z
% Returns:
%   d       - Divergence
dvx = px - circshift(px,1,1);  dvx(1,:,:) = px(1,:,:);  dvx(end,:,:) = -px(end-1,:,:);
dvy = py - circshift(py,1,2);  dvy(:,1,:) = py(:,1,:);  dvy(:,end,:) = -py(:,end-1,:);
dvz = pz - circshift(pz,1,3);  dvz(:,:,1) = pz(:,:,1);  dvz(:,:,end) = -pz(:,:,end-1);
d = dvx + dvy + dvz;

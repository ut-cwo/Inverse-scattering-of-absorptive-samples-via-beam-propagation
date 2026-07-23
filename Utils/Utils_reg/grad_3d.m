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
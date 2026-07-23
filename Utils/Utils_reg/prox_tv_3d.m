%% References
% A modified implementation of proximal 3D TV regularization as described 
% by Beck and Teboulle [1] adding a warm starting criterion and performing 
% minor adjustments for code performance with large gpuArrays.
%
% [1]   A. Beck and M. Teboulle. Fast gradient-based algorithms for 
%       constrained  total variation image denoising and deblurring 
%       problems. Image Processing, IEEE Transactions on, 
%       18(11):2419--2434, 2009
function [u, dual] = prox_tv_3d(x, gamma, numIters, warmStart, warmIts, dual)
%% prox_tv_3d: Proximal operator for 3D isotropic TV
% Arguments:
%   x       - Object of proximal TV
%   gamma   - Regularization strength
%   numIters- Number of tv iterations
%   dual    - Holds values for warm-start
% Returns:
%   u       - Regularized object
%   dual    - Write values for warm-start
arguments
    x           (:,:,:) {mustBeUnderlyingType(x,'single')}
    gamma
    numIters            = 100
    warmStart   logical = true
    warmIts             = 10
    dual        struct  = []
end

% Initialize
if ~isempty(dual)
    numIters    = warmIts;
    p = dual.p; q = dual.q; o = dual.o;
else
    p           = zeros(size(x), 'like', x); q = p; o = p;
end

% Dual step for 3D
c = 1/(12*gamma);

% Extrapolated/momentum points
r = p; s = q; k = o;

% Previous projected points
t = 1;

for it = 1:numIters
    % Operator evaluated at the momentum point
    [dx, dy, dz] = grad_3d(x - gamma * div_3d(r, s, k));

    % FGP momentum
    t_new   = (1 + sqrt(1 + 4*t^2))/2;
    beta    = (t - 1)/t_new;

    % Gradient-projection and momentum step
    [p, q, o, r, s, k] = arrayfun(@projMom,r,s,k,dx,dy,dz, ...
                                             c,p,q,o,beta);

    % Update
    t       = t_new;
end

% Final image reconstructed from the projected dual point
d       = div_3d(p, q, o);
u       = x - gamma * d;

if warmStart
    % Store computed dual variables for warm start
    dual(1).p  = p; dual(1).q = q; dual(1).o = o;
end

%% Array Function
function [pn, qn, on, rn, sn, kn] = projMom(r, s, k, dx, dy, dz, c, p, q, o, beta)
    % Gradient-projection step
    rt  = r - c*dx;
    st  = s - c*dy;
    kt  = k - c*dz;

    % Bound to isotropic L2 ball
    w   = max(sqrt(rt.^2 + st.^2 + kt.^2), 1);
    pn   = rt ./ w;
    qn   = st ./ w;
    on   = kt ./ w;

    % Momentum
    rn   = pn + beta*(pn - p);
    sn   = qn + beta*(qn - q);
    kn   = on + beta*(on - o);
end
end
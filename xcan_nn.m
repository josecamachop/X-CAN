function [P, P0, output] = xcan_nn(X,pcs,pen1,XtX,lambda1,XXt,lambda2,flag_nn, init)

R = max(pcs);

%% algorithm options
% if using lbfgsb (for nonnegativity constraints)
sz = size(X);
NN = sum(sz)*R + R + sz(2); % number of parameters (A, B, D, m)
l  = -inf(NN,1);    % lower bound
for n = 1:length(sz)
    if flag_nn(n)==1
        l(sum(sz(1:n-1))*R+1:sum(sz(1:n))*R) = zeros(sz(n)*R,1);
    end
end
l(end-R-sz(2)+1:end-sz(2)) = zeros(R,1); % constrain D to be non-negative
l(end-sz(2)+1:end) = zeros(sz(2),1); % constrain m to be non-negative
u  = inf(NN,1);              % there is no upper bound

options = struct('m', 5, 'printEvery', 100, 'maxIts', 10000,'maxTotalIts', 100000,'pgtol',1e-10,'lb',l,'ub',u,'flag_nn',flag_nn);

%% you need to give the above options as an input to xpca_opt; otherwise, it will just use the defaults.
[P, P0, output] = xcan_opt(X,pcs,pen1,lambda1, XtX, lambda2, XXt,'alg_options',options, 'alg','lbfgsb', 'init', init);
    
function [P, P0, output] = xcan_ncg(X,pcs,pen1,XtX,lambda1,XXt,lambda2, init)
% flag_nn is an array of size 2 indicating whether the factor matrices A
% and/or B are nonnegative. 

%% algorithm options
% if using ncg to fit the model 
options = ncg('defaults');
options.Display ='iter';
options.MaxIters      = 1e4;
options.MaxFuncEvals  = 1e5;
options.RelFuncTol    = 1e-10;
options.StopTol       = 1e-10;
options.DisplayIters  = 100;

%% you need to give the above options as an input to xpca_opt; otherwise, it will just use the defaults.
[P, P0, output] = xcan_opt(X,pcs,pen1,lambda1, XtX, lambda2, XXt,'alg_options',options, 'alg','ncg', 'init', init);


    
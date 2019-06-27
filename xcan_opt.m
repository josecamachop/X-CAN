function [P, P0, output] = xcan_opt(Z,pcs,pen1,lambda1, XtX, lambda2, XXt, varargin)
%SGPCA_OPT Fits SGPCA via optimization.
%
%   P = SGPCA_OPT(X,R) fits an R-component SGPCA model to Z.
%   The result P returns the factor matrices. lambda is the regularization parameter.
%   XX is the structural correlation, can be the correlation matrix, MEDA or similar.
%   The function being optimized is ....
%
%   K = SGPCA_OPT(X,R,'param',value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'alg' - Specfies optimization algorithm (default: 'ncg')
%      'ncg'   Nonlinear Conjugate Gradient Method
%      'lbfgs' Limited-Memory BFGS Method
%      'tn'    Truncated Newton
%
%   'init' - Initialization for factor matrices. (default:
%   'random'). This can be a cell array with the initial matrices or
%   one of the following strings:
%      'random' Randomly generated via randn function
%
%   'alg_options' - Parameter settings for selected optimization
%   algorithm. For example, type OPTIONS = NCG('defaults') to get
%   the NCG algorithm options which can then me modified as passed
%   through this function to NCG.
%
%   [P, P0] = SGPCA_OPT(...) also returns the initial guess.
%
%   [P, P0, OUT] = SGPCA_OPT(...) also returns a structure with the
%   optimization exit flag, the final relative fit, and the full
%   output from the optimization method. The fit is defined as 
%
%      FIT = 100 * (1 - ( F(K) / F(0) )).
%


%% Check for POBLANO
if ~exist('poblano_params','file')
    error(['SGPCA_OPT requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

%% Error checking
if (nargin < 6)
    error('Error: invalid input arguments');
end

%% Set parameters
params = inputParser;
params.addParamValue('alg', 'ncg', @(x) ismember(x,{'ncg','tn','lbfgs','lbfgsb'}));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addOptional('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        fhandle = @ncg;
    case 'tn'
        fhandle = @tn;
    case 'lbfgs'
        fhandle = @lbfgs;
    case 'lbfgsb'
        fhandle = @lbfgsb;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(fhandle, 'defaults');
else
    options = params.Results.alg_options;
end
        

%% Initialization (svd-based init can be added here!)
sz = size(Z);
N = length(sz);

P0=params.Results.init; 

withm = length(find(0==pcs));

%% Fit SGPCA using SGPCA_OPT
if strcmp(params.Results.alg,'lbfgsb')
    options.x0 = fac_to_vec(P0);
    [out.X, out.F, out.info_lbfgsb] = feval(fhandle, @(x)xcan_fun(x,Z,pen1, lambda1, XtX, lambda2, XXt, withm),options.lb, options.ub, options);
else 
    out  = feval(fhandle, @(x)xcan_fun(x,Z,pen1, lambda1, XtX, lambda2, XXt, withm), fac_to_vec(P0), options);
end
%normsqr =  norm(Z, 'fro')^2;

P = vec_to_fac(out.X, Z);
if nargout > 2
    if  strcmp(params.Results.alg,'lbfgsb')
         output.ExitFlag  = out.info_lbfgsb.lbfgs_message1;
         output.OptOut = out.info_lbfgsb;
         output.F = out.F;
    else
         output.ExitFlag  = out.ExitFlag;
         output.F = out.F;
    end
end



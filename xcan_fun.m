function [f,g] = xcan_fun(x, Z, pen1,lambda1, XtX, lambda2, XXt, withm)
%SGPCA_FUN Calculate function and gradient for SGPCA function.
%
%  [F,G] = SGPCA_FUN(x,Z, lambda, XX) where x is a vector containing the entries of the
%  components of the model, Z is the data matrix to be fit, lambda
%  regularization parameter and XX, e.g.,structural correlation, can be the
%  correlation matrix, MEDA or similar. 
%


%% Convert x to a cell array of matrices
A = vec_to_fac(x,Z);

%% Call cp_fit and cp_gradient using cp_fg
[f,g] = xcan_fg(Z, A, pen1, lambda1, XtX, lambda2, XXt, withm);




function XX = crossprod(X,opt)

% Cross-product matrix computation from data matrix X, following: 
% XX = X'*X./sqrt(diag(X'*X)*diag(X'*X)')
%
% XX = crossprod(X)    % minimum call
% XX = crossprod(X,opt)    % complete call
%
%
% INPUTS:
%
% X: [NxM] data set 
%
% opt: [1x1] options:
%       1: XX from X with no preprocessing, by default
%       2: XX from abs(X)
%       3: XX after baseline correction
%
%
% OUTPUTS:
%
% XX: [MxX] cross-product matrix.
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 25/Jun/19
%
% Copyright (C) 2019  University of Granada, Granada
% Copyright (C) 2019  Jose Camacho Paez
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 1, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
if nargin < 2 || isempty(opt), opt = 1; end;

%% Main code

switch opt,
    case 1, % No prep
        X = X; 
    case 2, % Absolute values
        X = abs(X);
    case 3, % Delete baseline
        X = X - ones(size(X,1),1)*min(X);
end

XX = X'*X;
XX = XX./sqrt(diag(XX)*diag(XX)');
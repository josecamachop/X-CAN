function [XP, XT, m] = xcan(Data,pcs,XtX,lc,XXt,lr,nn)
   
% XCAN algorithm
%
% [XP, XT] = xcan(Data,pcs)    % minimum call
% [XP, XT, m] = xcan(Data,pcs,XtX,lc,XXt,lr,nn)    % complete call
%
%
% INPUTS:
%
% Data: [NxM] data set 
%
% pcs: [1xA] Components considered (e.g. pcs = 1:2 selects the first two 
%   components). 
%
% XtX: [MxM] cross-product in the columns (ones by default).
%
% lc: [1x1] penalty coefficient in the columns (lc>=0, 0 by default).
%
% XXt: [NxN] cross-product in the rows (ones by default).
%
% lr: [1x1] penalty coefficient in the rows (ll>=0, 0 by default).
%
% nn: [2x1] non-negativity contraints in the rows and columns.
%   - [] or [0, 0]: contraints deactivated.
%   - [1 0]: non-neg in rows.
%   - [1 1]: non-neg in columns.
%   - [1 1]: non-neg in both rows and columns.
%
%
% OUTPUTS:
%
% XP: [MxA] loadings.
%
% XT: [NxA] scores.
%
% m: [1xM] baseline (only if pcs include the value 0).
%
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% coded by: Evrim Acar Ataman (evrim.acarataman@gmail.com)
% last modification: 26/Jun/19
%
% Copyright (C) 2019  University of Granada, Granada
% Copyright (C) 2019  Jose Camacho Paez
% Copyright (C) 2019  Evrim Acar Ataman 
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

% Correct data to Frob. norm 1 and mean center if 0 in pcs specified

%% Parameters checking

% Set default values
routine=dbstack;
assert (nargin >= 2, 'Error in the number of arguments. Type ''help %s'' for more info.', routine(1).name);
s=size(Data);
if nargin < 3 || isempty(XtX), XtX = ones(s(2)); end;
if nargin < 4 || isempty(lc), lc = 0; end;
if nargin < 5 || isempty(XXt), XXt = ones(s(1)); end;
if nargin < 6 || isempty(lr), lr = 0; end;
if nargin < 7 || isempty(nn), nn = []; end;


%% Main code

Nf = norm(Data,'fro');
Data = Data/Nf;

% Apply minimum threshold to XP matrices
thres1 = 1e-2;
ind = find(abs(XtX)<thres1); % Threshold for XtX
XtX(ind) = sign(XtX(ind)+1e-300*randn(size(XtX(ind))))*thres1;

thres2 = 1e-2;
ind = find(abs(XXt)<thres2); % Threshold for XXt
XXt(ind) = sign(XXt(ind)+1e-300*randn(size(XXt(ind))))*thres2;

% Build initial solution
P0 = {};
if find(pcs==0),
    m = mean(Data);
    Xc = Data - ones(size(Data,1),1)*m;
else
    m = zeros(1,s(2));
    Xc = Data;
end
P0{4}=m;

[P0{1},P0{3},P0{2}]=svd(Xc);
P0{1}=P0{1}(:,pcs(find(pcs>0)));
P0{2}=P0{2}(:,pcs(find(pcs>0)));
P0{3}=diag(P0{3}(pcs(find(pcs>0)),pcs(find(pcs>0))));
    
for j=pcs(find(pcs>0)),
    sPjind = find(abs(P0{1}(:,j))==max(abs(P0{1}(:,j))),1);
    sPj = sign(P0{1}(sPjind,j));
    P0{1}(:,j) = sPj*P0{1}(:,j);
    P0{2}(:,j) = sPj*P0{2}(:,j);
end
    
% Optimization
if isempty(nn) | isempty(find(nn>0)),
    P = xcan_ncg(Data,pcs,1,XtX,lc,XXt,lr,P0); 
else
    P = xcan_nn(Data,pcs,1,XtX,lc,XXt,lr,nn,P0);
end

% Computing outputs
D = P{3}*Nf;
U = P{1};
m = P{4}'*Nf;
XP = P{2};
XT = U*diag(D);

% Reorder according to corrected variance
Xc = Data*Nf - ones(size(Data,1),1)*m;
for j = pcs(find(pcs>0)),    
     Xp = (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*Xc*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
     varP(j) = trace(Xp'*Xp);
end

[kk,ord] = sort(varP,'Descend');
XP = XP(:,ord);
XT = XT(:,ord);

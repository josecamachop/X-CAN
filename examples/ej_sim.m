%% XCAN example of simulated data
% J. Camacho, E. Acar, M. Rasmussen, R. Bro. Cross-product Penalized 
% Component Analysis (XCAN), Submitted to Chemometrics and Intelligent 
% Laboratory Systems, 2019.
%
% Needs the MEDA Toolbox and the XCAN (path should be properly set)
%
% If you use these data, please add a reference to the above paper.

% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 27/Jun/19

close all
clear
clc

%% Simulated data

vars=5;
obs=5;

X=simuleMV(obs,vars,10);
X2=simuleMV(obs,vars,10);
X3=simuleMV(obs,vars,10);

X = [X zeros(5,5);zeros(5,5) X2;zeros(5,5) 2*X3] + 0.15*randn(15,10);

%save data/Simul X


%% Cross-product (XP) matrices

load data/Simul

XXt = crossprod(X');
plot_map(XXt)
ylabel('XXt','FontSize',20)

XX = crossprod(X);
plot_map(XX)
ylabel('XtX','FontSize',20)


%% XCAN: PCA vs XP-constraints

lambda = [0 .1 0 .1];
lambda2 = [0 0 .1 .1];
R = [3 3 3 3]
close all
clc

varX = trace(X'*X);
for i=1:length(lambda)
    
    [XP, XT]    = xcan(X,1:R(i),XX,lambda(i),XXt,lambda2(i));
    
    Xp = 0;
    for j=1:R(i),
        Xp = Xp + (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*X*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
    end
    varP = trace(Xp'*Xp);
    
    figure, subplot(2,1,1), title(sprintf('\\lambda_c: %g, \\lambda_r: %g, variance: %.0f%%',lambda(i),lambda2(i),100*varP/varX),'FontSize',20)
    
    hold on,  bar(XT), ylabel('Scores','FontSize',16), axis tight, subplot(2,1,2), bar(XP), ylabel('Loadings','FontSize',16), axis tight,
end


%% XCAN: Number of components

lambda = [.1 .1 .1];
lambda2 = [.1 .1 .1];
R = [1 2 3]
close all
clc

varX = trace(X'*X);
for i=1:length(lambda)
    
    [XP, XT]    = xcan(X,1:R(i),XX,lambda(i),XXt,lambda2(i));
        
        
    Xp = 0;
    for j=1:R(i),
        Xp = Xp + (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*X*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
    end
    varP = trace(Xp'*Xp);

    XP = [XP zeros(size(XP,1),3-i)];
    XT = [XT zeros(size(XT,1),3-i)];
    
    figure, subplot(2,1,1), title(sprintf('\\lambda_c: %g, \\lambda_r: %g, variance: %.0f%%',lambda(i),lambda2(i),100*varP/varX),'FontSize',20)
    
    hold on,  bar(XT), ylabel('Scores','FontSize',16), axis tight, subplot(2,1,2), bar(XP), ylabel('Loadings','FontSize',16), axis tight,
end


%% XCAN: Thresholded XP

load data/Simul

thr = 0.5; 
XXt = crossprod(X');
r = find(abs(XXt)<thr);
XXt(r) = 0;
plot_map(XXt)
ylabel('XXt','FontSize',20)

XX = crossprod(X);
r = find(abs(XX)<thr);
XX(r) = 0;
plot_map(XX)
ylabel('XtX','FontSize',20)

lambda = [.01];
lambda2 = [.01];
R = [3];
close all
clc

varX = trace(X'*X);
for i=1:length(lambda)
    
    [XP, XT]    = xcan(X,1:R(i),XX,lambda(i),XXt,lambda2(i));
    
    Xp = 0;
    for j=1:R(i),
        Xp = Xp + (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*X*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
    end
    varP = trace(Xp'*Xp);
    
    figure, subplot(2,1,1), title(sprintf('\\lambda_c: %g, \\lambda_r: %g, variance: %.0f%%',lambda(i),lambda2(i),100*varP/varX),'FontSize',20)
         
    hold on,  bar(XT), ylabel('Scores','FontSize',16), axis tight, subplot(2,1,2), bar(XP), ylabel('Loadings','FontSize',16), axis tight,  
end

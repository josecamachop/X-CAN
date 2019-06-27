%% XCAN example of olive data
% J. Camacho, E. Acar, M. Rasmussen, R. Bro. Cross-product Penalized 
% Component Analysis (XCAN), Submitted to Chemometrics and Intelligent 
% Laboratory Systems, 2019.
%
% Needs the MEDA Toolbox and the XCAN (path should be properly set)
%
% If you use these data, please add a reference to one of the below papers:
%
% de la Mata-Espinosa P, Bosque-Sendra JM, Bro R, Cuadros-Rodriguez L, 
% Discriminating olive and non-olive oils using HPLC-CAD and chemometrics, 
% Analytical and Bioanalytical Chemistry, 2011a, 399, 2083-2092.
%
% de la Mata-Espinosa P, Bosque-Sendra JM, Bro R, Cuadros-Rodriguez L, 
% Olive oil quantification of edible vegetable oil blends using 
% triacylglycerols chromatographic fingerprints and chemometric tools, 
% Talanta, 2011b, 85, 177-182.
%
% coded by: Jose Camacho Paez (josecamacho@ugr.es)
% last modification: 25/Jun/19

close all
clear
clc


%% Data viz

% load centered data, classes and wavelengths
load data/OliveData X class wavel

% PCA scores and loadings
scores_pca(X,1:2,[],0,[],[],class);
loadings_pca(X,1:2,0);

% load non-centered data and visualize it
load data/OliveData Xa

f=plot_vec2(Xa,{'Wavelength','Spectra (cm^{-1})'},class);
a=get(f,'Children');
set(a(1),'FontSize',22);
v = [1:13:length(Xa) length(Xa)];
set(a(1),'XTick',v);
set(a(1),'XTickLabel',wavel(v));
hold on, plot([27.5 27.5],[-0.01,0.03],'k--')
set(a(1),'XTickLabelRotation',45);


%% Cross-product (XP) matrices

% XP matrices after baseline correction
XXt = crossprod(Xa',3);
plot_map(XXt)
ylabel('XXt','FontSize',20)

XtX = crossprod(Xa,3);
plot_map(XtX)
ylabel('XtX','FontSize',20)

% XP matrix from classes
XXt2 = zeros(size(XXt));
for i=1:size(XXt,1),
    for j=1:size(XXt,1),
        if class(i) == class(j),
            XXt2(i,j) = 1;
        end
    end
end
plot_map(double(XXt2)); 
ylabel('XXt','FontSize',20)


%% XCAN: deactivated constraints Vs contraints in XtX

L=0.01;

lambdas = [0 L;0 0];
pcs = [6 6];

Data = Xa;
for i=1:size(lambdas,2)
    
    [XP, XT,m] = xcan(Data,0:pcs(i),XtX,lambdas(1,i),XXt,lambdas(2,i),[1,1]);
    
    plot_vec(m), ylabel('Loadings 0','FontSize',16),
    title('Baseline','FontSize',20)
    a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
    
    Data2 = Data - ones(size(Data,1),1)*m;
    varX = trace(Data2'*Data2);
    for j = 1:pcs(i),
        sPjind = find(abs(XT(:,j))==max(abs(XT(:,j))),1);
        sPj = sign(XT(sPjind,j));
        
        Xp = (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*Data2*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
        varP = trace(Xp'*Xp);
        
        plot_vec(sPj*XT(:,j),[],class);
        ylabel('Scores','FontSize',16),
        title(sprintf('XC: %i, variance: %.2f%%',j,100*(varP)/varX),'FontSize',20)
        a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
        
        plot_vec(sPj*XP(:,j)), ylabel('Loadings','FontSize',16),
        a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
    end
   
    E = Data - ones(size(Data,1),1)*m - XT(:,1:4)*XP(:,1:4)';        
    varE = trace(E'*E);
    ResidualVar = 100*(varE)/varX
end

        
%% From the classes

lambdas = [L;L];
pcs = 6;

Data = Xa;
for i=1:size(lambdas,2)
    
    [XP, XT,m] = xcan(Data,0:pcs(i),XtX,lambdas(1,i),XXt2,lambdas(2,i),[1,1]);
       
    plot_vec(m), ylabel('Loadings 0','FontSize',16),
    title('Baseline','FontSize',20)
    a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
    
    Data2 = Data - ones(size(Data,1),1)*m;
    varX = trace(Data2'*Data2);
    for j = 1:pcs(i),
        sPjind = find(abs(XT(:,j))==max(abs(XT(:,j))),1);
        sPj = sign(XT(sPjind,j));
        
        Xp = (XT(:,j)*pinv(XT(:,j)'*XT(:,j))*XT(:,j)')*Data2*(XP(:,j)*pinv(XP(:,j)'*XP(:,j))*XP(:,j)');
        varP = trace(Xp'*Xp);

        plot_vec(sPj*XT(:,j),[],class);
        ylabel('Scores','FontSize',16),
        title(sprintf('XC: %i, variance: %.2f%%',j,100*(varP)/varX),'FontSize',20)
        a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
       
        plot_vec(sPj*XP(:,j)), ylabel('Loadings','FontSize',16),
        a=get(gcf,'Children'); set(a,'PlotBoxAspectRatio',[1 .3 .3])
        
    end
    
    E = Data - ones(size(Data,1),1)*m - XT(:,1:4)*XP(:,1:4)';        
    varE = trace(E'*E);
    ResidualVar = 100*(varE)/varX
end



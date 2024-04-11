%% Uncoupling Preprocessing Workflow
% Load in your datasets
% FDG Data
load("LOAD2_HFD_FDG.mat")
% PTSM Data
load("LOAD2_HFD_PTSM.mat")

%% Male Data 
% FDG
% 18 month LOAD2 males
A = LOAD2_HFD_FDG{4,3};
% 18 HFD month LOAD2 males
B = LOAD2_HFD_FDG{7,3};

uncoupling_preprocessing(A,B)
load("preprocessed.mat")

M_FDG_ref = z_ref; M_FDG_ref_mean = z_ref_mean; M_FDG_exp2ref = z_exp2ref; M_FDG_exp2ref_mean = z_exp2ref_mean; M_FDG_exp2ref_SEM = z_exp2ref_SEM;

writetable(t,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Male_18mo_CDvHFD_FDG")
[rw,~] = size(t(:,2));
for i=1:rw
    M_FDG_sharedz_p05(i,1) = t(i,2); %#ok<*SAGROW>
end
M_FDG_sharedz_p05 = table2array(M_FDG_sharedz_p05);

% PTSM 
A2 = LOAD2_HFD_PTSM{4,3};
B2 = LOAD2_HFD_PTSM{7,3};
uncoupling_preprocessing(A2,B2)
load("preprocessed.mat")
M_PTSM_ref = z_ref; M_PTSM_ref_mean = z_ref_mean; M_PTSM_exp2ref = z_exp2ref; M_PTSM_exp2ref_mean = z_exp2ref_mean; M_PTSM_exp2ref_SEM = z_exp2ref_SEM;

writetable(t,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Male_18mo_CDvHFD_PTSM")
[rw,~] = size(t(:,2));
for i=1:rw
    M_PTSM_sharedz_p05(i,1) = t(i,2);
end
M_PTSM_sharedz_p05 = table2array(M_PTSM_sharedz_p05);
for i=1:rw
    [~,M_Uncoupled(i,1)] = ttest2(M_FDG_exp2ref(i,:),M_PTSM_exp2ref(i,:));
    if M_Uncoupled(i,1)<0.05
        M_Uncoupled_p05(i,1) = M_Uncoupled(i,1);
    else
        M_Uncoupled_p05(i,1) = 0;
    end
end

%% Female Data
% FDG
C = LOAD2_HFD_FDG{4,2};
D = LOAD2_HFD_FDG{7,2};
uncoupling_preprocessing(C,D)
load("preprocessed.mat")
F_FDG_ref = z_ref; F_FDG_ref_mean = z_ref_mean; F_FDG_exp2ref = z_exp2ref; F_FDG_exp2ref_mean = z_exp2ref_mean; F_FDG_exp2ref_SEM = z_exp2ref_SEM;

writetable(t,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Female_18mo_CDvHFD_FDG")

[rw,~] = size(t(:,2));
for i=1:rw
    F_FDG_sharedz_p05(i,1) = t(i,2);
end

F_FDG_sharedz_p05 = table2array(F_FDG_sharedz_p05);

% PTSM 
C2 = LOAD2_HFD_PTSM{4,2};
D2 = LOAD2_HFD_PTSM{7,2};
uncoupling_preprocessing(C2,D2)
load("preprocessed.mat")
F_PTSM_ref = z_ref; F_PTSM_ref_mean = z_ref_mean; F_PTSM_exp2ref = z_exp2ref; F_PTSM_exp2ref_mean = z_exp2ref_mean; F_PTSM_exp2ref_SEM = z_exp2ref_SEM;
writetable(t,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Female_18mo_CDvHFD_PTSM")

[rw,~] = size(t(:,2));
for i=1:rw
    F_PTSM_sharedz_p05(i,1) = t(i,2);
end

F_PTSM_sharedz_p05 = table2array(F_PTSM_sharedz_p05);

for i=1:rw
    [~,F_Uncoupled(i,1)] = ttest2(F_FDG_exp2ref(i,:),F_PTSM_exp2ref(i,:));
    if F_Uncoupled(i,1)<0.05
        F_Uncoupled_p05(i,1) = F_Uncoupled(i,1);
    else
        F_Uncoupled_p05(i,1) = 0;
    end
end

%% Z-Scores
% Take mean z-scores for comparator groups and write them to an Excel table
% Comparison of global z-scoring between reference and experimental cohorts on a
% per-region basis
zscores_sharedz = [M_FDG_ref_mean; M_FDG_exp2ref_mean; F_FDG_ref_mean; F_FDG_exp2ref_mean; M_PTSM_ref_mean; M_PTSM_exp2ref_mean; F_PTSM_ref_mean;F_PTSM_exp2ref_mean]; zscores_sharedz = transpose(zscores_sharedz);
z = array2table(zscores_sharedz,"RowNames",{'AI','AuDMV','CPu','Cg','CC','DLO','DLIVEnt','DI','ECT','Fornix','FrA','HIP','LO','MO','PtPR','PtA','PRH','PrL','M1','S1','RSC','M2','S2','TeA','TH','VO','V1V2'}, ...
    "VariableNames",{'M_FDG_REF','M_FDG_EXP2REF','F_FDG_REF','F_FDG_EXP2REF','M_PTSM_REF','M_PTSM_EXP2REF','F_PTSM_REF','F_PTSM_EXP2REF'});

writetable(z,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","ZScores")

%% Thresholded p-values for atlas modeler
pvals_sharedz = [M_FDG_sharedz_p05 M_PTSM_sharedz_p05 M_Uncoupled_p05 F_FDG_sharedz_p05 F_PTSM_sharedz_p05 F_Uncoupled_p05];
[rw,cl] = size(pvals_sharedz);
for i=1:rw
    for j=1:cl
        if pvals_sharedz(i,j)<0.05
            pvals_sharedz_thr(i,j) = pvals_sharedz(i,j);
        else
            pvals_sharedz_thr(i,j) = 0;
        end
    end
end
clear rw cl

pvals_sharedz = [M_FDG_sharedz_p05 M_PTSM_sharedz_p05 M_Uncoupled_p05 F_FDG_sharedz_p05 F_PTSM_sharedz_p05 F_Uncoupled_p05];
[rw,cl] = size(pvals_sharedz);
for i=1:rw
    for j=1:cl
        if pvals_sharedz(i,j)<0.05
            pvals_sharedz_thr(i,j) = pvals_sharedz(i,j);
        else
            pvals_sharedz_thr(i,j) = 0;
        end
    end
end
clear rw cl

z2 = array2table(pvals_sharedz_thr,"RowNames",{'AI','AuDMV','CPu','Cg','CC','DLO','DLIVEnt','DI','ECT','Fornix','FrA','HIP','LO','MO','PtPR','PtA','PRH','PrL','M1','S1','RSC','M2','S2','TeA','TH','VO','V1V2'}, ...
    "VariableNames",{'M_FDG','M_PTSM','M_Uncoupled','F_FDG','F_PTSM','F_Uncoupled'});

writetable(z2,"Uncoupling_18moHFD_REF_18moCD.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","P-vals<0.05")

%% Visualize Uncoupling Data

% Determine max, min values for setting axis limits
male_mean_range = [M_PTSM_exp2ref_mean M_FDG_exp2ref_mean];
female_mean_range = [F_PTSM_exp2ref_mean F_FDG_exp2ref_mean];
combined_mean_range = [male_mean_range female_mean_range];

min_m = min(male_mean_range); max_m = max(male_mean_range);
if abs(min_m)>abs(max_m)
    axis_lim_m = abs(min_m-.5);
else
    axis_lim_m = abs(max_m+.5);
end
min_f = min(female_mean_range); max_f = max(female_mean_range);
if abs(min_f)>abs(max_f)
    axis_lim_f = abs(min_f-.5);
else
    axis_lim_f = abs(max_f+.5);
end
abs_min = min(combined_mean_range); abs_max = max(combined_mean_range);
if abs(abs_min)>abs(abs_max)
    axis_lim_abs = abs(abs_min-.5);
else
    axis_lim_abs = abs(abs_max+.5);
end
% Generate Scatter Plots of Male Uncoupling, Female Uncoupling, and both
% Male
f1 = figure('units','inches','PaperPositionMode','auto','Name',"Male Uncoupling");
d = scatter(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,'filled','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[.3 .3 1]); distfromzero = sqrt(M_PTSM_exp2ref_mean.^2 + M_FDG_exp2ref_mean.^2); d.AlphaData = distfromzero; 
d.MarkerFaceAlpha = 'flat'; hold on
eb(1) = errorbar(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,M_PTSM_exp2ref_SEM, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,M_FDG_exp2ref_SEM, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [.5 .5 .5], 'LineWidth', .5)
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-axis_lim_m axis_lim_m -axis_lim_m axis_lim_m])
Xlm = xlim;
Ylm = ylim;
Xlb = mean(Xlm);
xlabel("X-Axis: Cerebral Perfusion Z-Score Relative to Control","Position",[Xlb (Ylm(1)-.5)],"VerticalAlignment","bottom","HorizontalAlignment","center","FontSize",8);
ylabel("Y-Axis: Cerebral Metabolic Uptake Z-Score Relative to Control","Position",[Xlm(1) mean(Ylm)],"VerticalAlignment","bottom","HorizontalAlignment","center","Rotation",90,"FontSize",8)
legend(d,{'Male'},"Location","northwest")
legend('boxoff')

% Female
f2 = figure('units','inches','PaperPositionMode','auto','Name',"Female Uncoupling");
d = scatter(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,'filled','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 .3 .3]);
distfromzero = sqrt(F_PTSM_exp2ref_mean.^2 + F_FDG_exp2ref_mean.^2); d.AlphaData = distfromzero; d.MarkerFaceAlpha = 'flat'; hold on
eb(1) = errorbar(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,F_PTSM_exp2ref_SEM, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,F_FDG_exp2ref_SEM, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [.5 .5 .5], 'LineWidth', .5)
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-axis_lim_f axis_lim_f -axis_lim_f axis_lim_f])
Xlm = xlim;
Ylm = ylim;
Xlb = mean(Xlm);
xlabel("X-Axis: Cerebral Perfusion Z-Score Relative to Control","Position",[Xlb (Ylm(1)-.5)],"VerticalAlignment","bottom","HorizontalAlignment","center","FontSize",8);
ylabel("Y-Axis: Cerebral Metabolic Uptake Z-Score Relative to Control","Position",[Xlm(1) mean(Ylm)],"VerticalAlignment","bottom","HorizontalAlignment","center","Rotation",90,"FontSize",8)
legend(d,{'Female'},"Location","northwest")
legend('boxoff')


% Both
f3 = figure(Units='inches',PaperPositionMode='auto',Name="Both Sexes Uncoupling");
d = scatter(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,'filled','MarkerEdgeColor',[0 0 1],'MarkerFaceColor',[.3 .3 1]);
distfromzero = sqrt(M_PTSM_exp2ref_mean.^2 + M_FDG_exp2ref_mean.^2); d.AlphaData = distfromzero; d.MarkerFaceAlpha = 'flat'; hold on
eb(1) = errorbar(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,M_PTSM_exp2ref_SEM, 'horizontal', 'LineStyle', 'none');
eb(2) = errorbar(M_PTSM_exp2ref_mean,M_FDG_exp2ref_mean,M_FDG_exp2ref_SEM, 'vertical', 'LineStyle', 'none');
set(eb, 'color', [.5 .5 .5], 'LineWidth', .5); hold on;
e = scatter(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,'filled','MarkerEdgeColor',[1 0 0],'MarkerFaceColor',[1 .3 .3]);
distfromzero = sqrt(F_PTSM_exp2ref_mean.^2 + F_FDG_exp2ref_mean.^2); e.AlphaData = distfromzero; e.MarkerFaceAlpha = 'flat'; hold on
eb2(1) = errorbar(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,F_PTSM_exp2ref_SEM, 'horizontal', 'LineStyle', 'none');
eb2(2) = errorbar(F_PTSM_exp2ref_mean,F_FDG_exp2ref_mean,F_FDG_exp2ref_SEM, 'vertical', 'LineStyle', 'none');
set(eb2, 'color', [.5 .5 .5], 'LineWidth', .5)
set(gca, 'XAxisLocation', 'origin', 'YAxisLocation', 'origin')
axis([-axis_lim_abs axis_lim_abs -axis_lim_abs axis_lim_abs])
Xlm = xlim;
Ylm = ylim;
Xlb = mean(Xlm);
Ylb = .99*Ylm(1);
xlabel("X-Axis: Cerebral Perfusion Z-Score Relative to Control","Position",[Xlb (Ylm(1)-.5)],"VerticalAlignment","bottom","HorizontalAlignment","center","FontSize",8);
ylabel("Y-Axis: Cerebral Metabolic Uptake Z-Score Relative to Control","Position",[Xlm(1) mean(Ylm)],"VerticalAlignment","bottom","HorizontalAlignment","center","Rotation",90,"FontSize",8)
legend([d e],{'Male','Female'},"Location","northwest")
legend('boxoff')

% Save outputs
exportgraphics(f1,[pwd '\Male_Output.png'],"BackgroundColor","white","Resolution",600)
exportgraphics(f2,[pwd '\Female_Output.png'],"BackgroundColor","white","Resolution",600)
exportgraphics(f3,[pwd '\MF_Output.png'],"BackgroundColor","white","Resolution",600)

%% Run connectomics on FDG and PTSM (Shared Distribution)
addpath(genpath(pwd))
% Load V2 which has groupings of nodes into systems/networks
load('roiLabels_wGroupings.mat')

% This will load colormaps contained in the package
load('colormaps.mat')
pval = 0.05;
zdata2{1,1} = 'LOAD2 18m HFD Zscores Shared Dist'; zdata2{2,1} = 'FDG ref'; zdata2{3,1} = 'PTSM ref'; zdata2{4,1} = 'FDG exp2ref'; zdata2{5,1} = 'PTSM exp2ref'; zdata2{1,2} = 'Male'; zdata2{1,3} = 'Female';
zdata2{2,2} = M_FDG_ref; zdata2{2,3} = F_FDG_ref; zdata2{3,2} = M_PTSM_ref; zdata2{3,3} = F_PTSM_ref; zdata2{4,2} = M_FDG_exp2ref; zdata2{4,3} = F_FDG_exp2ref; zdata2{5,2} = M_PTSM_exp2ref; zdata2{5,3} = F_PTSM_exp2ref;
save("Uncoupling_ZData_18Mo.mat","zdata2")
covariance_analysis_tier1(zdata2,roi_labels,bluered_cmap,pval)

%% Run isomorphic analyses using Jaccard distance as metric to see the coupling differences between PTSM and FDG
% Load Tier1 Ouput of FDG and PTSM Shared Distribution Z-Score Data

% Compute Jaccard Distance Between Graphs
% JaccardPET(

%% Run connectomics on FDG and PTSM Z-Scored Data (Independent Distributions)
% Create new nested cell array of z-scored FDG and PTSM Data
zdata{1,1} = 'LOAD2 18m HFD Zscores'; zdata{2,1} = 'FDG ref'; zdata{3,1} = 'PTSM ref'; zdata{4,1} = 'FDG exp'; zdata{5,1} = 'PTSM exp'; zdata{1,2} = 'Male'; zdata{1,3} = 'Female';
zdata{2,2} = M_FDG_ref; zdata{2,3} = F_FDG_ref; zdata{3,2} = M_PTSM_ref; zdata{3,3} = F_PTSM_ref; zdata{4,2} = M_FDG_exp; zdata{4,3} = F_FDG_exp; zdata{5,2} = M_PTSM_exp; zdata{5,3} = F_PTSM_exp;

% Prepare the necessary inputs
% Add code package to matlab path. (THIS PATH WILL VARY FOR EACH USER)
addpath(genpath(pwd))
% Load V2 which has groupings of nodes into systems/networks
load('roiLabels_wGroupings.mat')

% This will load colormaps contained in the package
load('colormaps.mat')
pval = 0.05;
%covariance_analysis_tier1(zdata,roi_labels,bluered_cmap,pval)

%% Compute the Jaccard Index as a similarity metric between different uncoupling covariance matrices



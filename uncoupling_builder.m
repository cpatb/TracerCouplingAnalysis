function uncoupling_builder(FDG,PTSM)
% Takes in covariance matrix arrays from outputs of tier1 connectomics
% analyses of FDG and PTSM data, respectively
narginchk(2,3)
[rwF, clF] = size(FDG);
[rwP, clP] = size(PTSM);
% cell structure check variable names in first row col
 FDG_rNames = FDG(:,1);    % row subgroup labels FDG
 FDG_cNames = FDG(1,:);    % column subgroup labels FDG
 PTSM_rNames = PTSM(:,1); % row subgroup labels PTSM
 PTSM_cNames = PTSM(1,:); % column subgroup labels PTSM

 % there is a +1 for (1,1) which is why most loops start at 2.
 NrF = length(FDG_rNames);       % number of row groupings FDG
 NcF = length(FDG_cNames);       % number of column groupings FDG
 NrP = length(PTSM_rNames);       % number of row groupings PTSM
 NcP = length(PTSM_cNames);       % number of column groupings PTSM

if rwF ~= rwP
    fprintf(2,'Number of rows varies between tracers. Check that input cell arrays were made correctly.\n')
    return % exit
elseif clF ~= clP
    fprintf(2,'Number of columns varies between tracers. Check that input cell arrays were made correctly.\n')
    return
end

%% Create sex*treatment supermatrices
for i=2:rwF
    for j=2:clF
        FDG_data = FDG{i,j};
        PTSM_data = PTSM{i,j};
        [region_countF,cohort_countF] = size(FDG_data);
        [region_countP,cohort_countP] = size(PTSM_data);
        % run mean imputation
        if region_countF == region_countP
            if cohort_countF > cohort_countP
                for k=1:region_countP
                    for l=cohort_countP:cohort_countF
                        PTSM_data(k,l) = mean(PTSM_data(k,:));
                    end
                end
                %PTSM_data(:,end+1:cohort_countF)=missing;
                %super = NaN(2*region_countF,cohort_countF);
            elseif cohort_countP > cohort_countF
                for m=1:region_countF
                    for n=cohort_countF:cohort_countP
                        FDG_data(m,n) = mean(FDG_data(m,:));
                    end
                end
                FDG_data(:,end+1:cohort_countP)=missing;
                %super = NaN(2*region_countP,cohort_countP);
            end
        else
            fprintf(2,'Number of regions varies between tracers. Check that nested matrices were made correctly.\n')
        return
        end
        super{i,j} = [FDG_data; PTSM_data];
        supert = transpose(super{i,j});
        % Initialize the covariance supermatrix
        covariance_supermatrix = zeros(2*region_countF, 2*region_countF);

        % Compute covariances for FDG-PTSM interactions
        for o = 1:(region_countF*2)
            for p = 1:(region_countF*2)
                covariance_supermatrix(o,p) = corr(supert(:,o), supert(:,p),"type","Pearson");
            end
        end
        uncoupling_supermatrix{i,j} = covariance_supermatrix;
        clear covariance_supermatrix
        uncoupling_supermatrix{1,j} = FDG{1,j};
        
    end
    uncoupling_supermatrix{i,1} = FDG{i,1}; 
end
uncoupling_supermatrix{1,1} = "FDG_PTSM_Supermatrices";

%% Isolate, transform, run connectomics on uncoupling submatrix
% Isolate uncoupling matrix (quadrants I and III)

% Sizes of the matrices
rowlim_FDG_FDG = 27;
collim_FDG_FDG = 27;
rowlim_PTSM_PTSM = 54;
collim_PTSM_PTSM = 54;
rowlim_FDG_PTSM = 54;
collim_FDG_PTSM = 27;

% Extract the submatrices using specific sizes
[rw_super, cl_super] = size(uncoupling_supermatrix);
for i=2:rw_super
    for j=2:cl_super
        fullsuper = uncoupling_supermatrix{i,j}; FDG_FDG = fullsuper(1:rowlim_FDG_FDG, 1:collim_FDG_FDG);
        FDG_PTSM = fullsuper(rowlim_FDG_FDG+1:rowlim_FDG_PTSM, 1:collim_FDG_PTSM);
        PTSM_PTSM = fullsuper(rowlim_FDG_FDG+1:rowlim_PTSM_PTSM, collim_FDG_FDG+1:collim_PTSM_PTSM);
        FDG_FDG_matrix{i,j} = FDG_FDG; FDG_PTSM_matrix{i,j} = FDG_PTSM; PTSM_PTSM_matrix{i,j} = PTSM_PTSM;
        FDG_FDG_matrix(1,j) = uncoupling_supermatrix(1,j);FDG_PTSM_matrix(1,j) = uncoupling_supermatrix(1,j); PTSM_PTSM_matrix(1,j) = uncoupling_supermatrix(1,j);
    end
    FDG_FDG_matrix(i,1) = uncoupling_supermatrix(i,1);
    FDG_PTSM_matrix(i,1) = uncoupling_supermatrix(i,1); FDG_PTSM_matrix{1,1} = "Coupled";
    PTSM_PTSM_matrix(i,1) = uncoupling_supermatrix(i,1);
end

% Transform uncoupling matrix
[rw_uncpld,cl_uncpld] = size(FDG_PTSM_matrix);
for i=2:rw_uncpld
    for j=2:cl_uncpld
        transformed_uncoupled{i,j} = FDG_PTSM_matrix{i,j}*-1;
        transformed_uncoupled{1,j} = FDG_PTSM_matrix{1,j};
    end
    transformed_uncoupled{i,1} = FDG_PTSM_matrix{i,1}; transformed_uncoupled{1,1} = "Transformed + Coupled";
end
%% Run Tier1 Analysis on Uncoupling Matrices
% Permutation testing of significant edge covariance 
% 10,000 permutations
% We use FDG structure, since we have already proven that FDG and PTSM
% datasets have the same structure
cov_permP = FDG_rNames;    cov_permP(1,1:NcF) = FDG_cNames;    % permutation p value matrices
for r=2:NrF % every row
    for c=2:NcF % every column
        [~,~,cov_permP{r,c}] = randshiftnull_cov(super{r,c},10000);
    end
end
% Extract the permutation uncoupling submatrix using specific sizes
[rw_perm, cl_perm] = size(cov_permP);
for i=2:rw_perm
    for j=2:cl_perm
        fullsuper = cov_permP{i,j}; 
        FDG_PTSM = fullsuper(rowlim_FDG_FDG+1:rowlim_FDG_PTSM, 1:collim_FDG_PTSM);
        FDG_PTSM_perm_matrix{i,j} = FDG_PTSM;
        FDG_PTSM_perm_matrix(1,j) = uncoupling_supermatrix(1,j);
    end
    FDG_PTSM_perm_matrix(i,1) = uncoupling_supermatrix(i,1); FDG_PTSM_perm_matrix{1,1} = "Coupled";
end
% Create thresholded covariance matrices
cnt=0;
for r=2:Nr
    for c=2:Nc
        cnt=cnt+1;
        mask = logical(FDG_PTSM_perm_matrix{r,c}>0.05);     % binary mask of NONsignificant edges
        mat = transformed_uncoupled{r,c};                          % pull out covariance matrix
        mat(mask)=0;                                % zero edges greater than pval
        mat(logical(eye(length(mat))))=0;           % zero the diagonal (self connections)
        thr_cov{r,c} = mat;                   % store thresholded covariance matrix
        clear mask mat
    end
end



% %% Run connectomics on transformed uncoupling matrices
% load('roiLabels_wGroupings.mat','roi_labels')
% load('colormaps.mat','bluered_cmap')
% N = length(roi_labels);
% for i=2:rw_transformed
%     for j=2:cl_transformed
%         matrix = transformed_uncoupled{i,j};
%         [agrMat, mrccPartition,allPartitions] = mrcc_wrapper(matrix*-1,10000,1,roi_labels);
%         mrccPartition = fcn_relabel_partitions(mrccPartition);
%         [X,Y,idx_ord]=grid_communities(mrccPartition);
%         f(i,j)= figure('units','inches','position',[1 1 6 6],'paperpositionmode','auto');
%         imagesc(matrix(idx_ord,idx_ord)); axis square
%         colormap(bluered_cmap); clim([-1 1])
%         hold on
%         plot(X,Y,'k','LineWidth',2)
%         xticks(1:1:N); yticks(1:1:N);
%         yticklabels(roi_labels(idx_ord)); xticklabels(roi_labels(idx_ord)); xtickangle(90)
%         colorbar
%         title(strcat(transformed_uncoupled{i,1}," ",transformed_uncoupled{1,j}," Community Ordered Covariance"))
%         exportgraphics(f(i,j),strcat(pwd,"\",transformed_uncoupled{i,1},"_",transformed_uncoupled{1,j},"_coupling.png"),"Resolution",300)
%     end
% end
% close all



save("output.mat","uncoupling_supermatrix","transformed_uncoupled","FDG_PTSM_matrix")
close all
end
function uncoupling_builder(T1,T2)
% Takes in covariance matrix arrays from outputs of tier1 connectomics
% analyses of Tracer 1 (T1) and Tracer 2 (T2 )data, respectively and returns connectomic analyses of a 2-tracer coupled covariance matrix
narginchk(2,2)
[rT1, cT1] = size(T1);
[rT2, cT2] = size(T2);
% cell structure check variable names in first row col
 T1_rNames = T1(:,1);    % row subgroup labels T1
 T1_cNames = T1(1,:);    % column subgroup labels T1
 %T2_rNames = T2(:,1); % row subgroup labels T2
 %T2_cNames = T2(1,:); % column subgroup labels T2

 % there is a +1 for (1,1) which is why most loops start at 2.
 NrT1 = length(T1_rNames);       % number of row groupings T1
 NcT1 = length(T1_cNames);       % number of column groupings T1
 %NrT2 = length(T2_rNames);       % number of row groupings T2
 %NcT2 = length(T2_cNames);       % number of column groupings T2

if rT1 ~= rT2
    fprintf(2,'Number of rows varies between tracers. Check that input cell arrays were made correctly.\n')
    return % exit
elseif cT1 ~= cT2
    fprintf(2,'Number of columns varies between tracers. Check that input cell arrays were made correctly.\n')
    return
end

%% Create sex*treatment supermatrices
for i=2:rT1
    for j=2:clF
        T1_data = T1{i,j};
        T2_data = T2{i,j};
        [region_ct_T1,cohort_ct_T1] = size(T1_data);
        [region_ct_T2,cohort_ct_T2] = size(T2_data);
        % run mean imputation - - - - - - NOTE THAT THIS IS FOR FIRST PASS ESTIMATIONS ONLY
        if region_ct_T1 == region_ct_T2
            if cohort_ct_T1 > cohort_ct_T2
                for k=1:region_ct_T2
                    for l=cohort_ct_T2:cohort_ct_T1
                        T2_data(k,l) = mean(T2_data(k,:));
                    end
                end
            elseif cohort_ct_T2 > cohort_ct_T1
                for m=1:region_ct_T1
                    for n=cohort_ct_T1:cohort_ct_T2
                        T1_data(m,n) = mean(T1_data(m,:));
                    end
                end
                T1_data(:,end+1:cohort_ct_T2)=missing;
            end
        else
            fprintf(2,'Number of regions varies between tracers. Check that nested matrices were made correctly.\n')
        return
        end
        super{i,j} = [T1_data; T2_data]; %#ok<AGROW>
        supert = transpose(super{i,j});
        % Initialize the covariance supermatrix
        covariance_supermatrix = zeros(2*region_ct_T1, 2*region_ct_T1);

        % Compute covariances for T1-T2 interactions
        for o = 1:(region_ct_T1*2)
            for p = 1:(region_ct_T1*2)
                covariance_supermatrix(o,p) = corr(supert(:,o), supert(:,p),"type","Pearson");
            end
        end
        uncoupling_supermatrix{i,j} = covariance_supermatrix; %#ok<AGROW>
        clear covariance_supermatrix
        uncoupling_supermatrix{1,j} = T1{1,j}; %#ok<AGROW>
        
    end
    uncoupling_supermatrix{i,1} = T1{i,1}; 
end
uncoupling_supermatrix{1,1} = "T1_T2_Supermatrices";

%% Isolate, transform, run connectomics on uncoupling submatrix
% Isolate uncoupling matrix (quadrants I and III)

% Sizes of the matrices
rowlim_T1_T1 = rT1;
collim_T1_T1 = rT1;
rowlim_T2_T2 = rT1+rT2;
collim_T2_T2 = rT1+rT2;
rowlim_T1_T2 = rowlim_T2_T2;
collim_T1_T2 = collim_T1_T1;

% Extract the submatrices using specific sizes
[rw_super, cl_super] = size(uncoupling_supermatrix);
for i=2:rw_super
    for j=2:cl_super
        fullsuper = uncoupling_supermatrix{i,j}; T1_T1 = fullsuper(1:rowlim_T1_T1, 1:collim_T1_T1);
        T1_T2 = fullsuper(rowlim_T1_T1+1:rowlim_T1_T2, 1:collim_T1_T2);
        T2_T2 = fullsuper(rowlim_T1_T1+1:rowlim_T2_T2, collim_T1_T1+1:collim_T2_T2);
        T1_T1_matrix{i,j} = T1_T1; T1_T2_matrix{i,j} = T1_T2; T2_T2_matrix{i,j} = T2_T2; %#ok<AGROW>
        T1_T1_matrix(1,j) = uncoupling_supermatrix(1,j);T1_T2_matrix(1,j) = uncoupling_supermatrix(1,j); T2_T2_matrix(1,j) = uncoupling_supermatrix(1,j); %#ok<AGROW>
    end
    T1_T1_matrix(i,1) = uncoupling_supermatrix(i,1); T1_T1_matrix{1,1} = T1{1,1};
    T1_T2_matrix(i,1) = uncoupling_supermatrix(i,1); T1_T2_matrix{1,1} = "Coupled";
    T2_T2_matrix(i,1) = uncoupling_supermatrix(i,1); T2_T2_matrix{1,1} = T2{1,1};
end

% Transform uncoupling matrix
[rw_uNcT2ld,cl_uNcT2ld] = size(T1_T2_matrix);
for i=2:rw_uNcT2ld
    for j=2:cl_uNcT2ld
        transformed_uncoupled{i,j} = T1_T2_matrix{i,j}*-1; %#ok<AGROW>
        transformed_uncoupled{1,j} = T1_T2_matrix{1,j}; %#ok<AGROW>
    end
    transformed_uncoupled{i,1} = T1_T2_matrix{i,1}; transformed_uncoupled{1,1} = "Transformed + Coupled";
end
%% Run Tier1 Analysis on Uncoupling Matrices
% Permutation testing of significant edge covariance 
% 10,000 permutations
% We use T1 structure, since we have already proven that T1 and T2
% datasets have the same structure
cov_permP = T1_rNames;    cov_permP(1,1:NcT1) = T1_cNames;    % permutation p value matrices
for r=2:NrT1 % every row
    for c=2:NcT1 % every column
        [~,~,cov_permP{r,c}] = randshiftnull_cov(super{r,c},10000);
    end
end
% Extract the permutation uncoupling submatrix using specific sizes
[rw_perm, cl_perm] = size(cov_permP);
for i=2:rw_perm
    for j=2:cl_perm
        fullsuper = cov_permP{i,j}; 
        T1_T2 = fullsuper(rowlim_T1_T1+1:rowlim_T1_T2, 1:collim_T1_T2);
        T1_T2_perm_matrix{i,j} = T1_T2;
        T1_T2_perm_matrix(1,j) = uncoupling_supermatrix(1,j);
    end
    T1_T2_perm_matrix(i,1) = uncoupling_supermatrix(i,1); T1_T2_perm_matrix{1,1} = "Coupled";
end
% Create thresholded covariance matrices
cnt=0;
for r=2:Nr
    for c=2:Nc
        cnt=cnt+1;
        mask = logical(T1_T2_perm_matrix{r,c}>0.05);     % binary mask of NONsignificant edges
        mat = transformed_uncoupled{r,c};                          % pull out covariance matrix
        mat(mask)=0;                                % zero edges greater than pval
        mat(logical(eye(length(mat))))=0;           % zero the diagonal (self connections)
        thr_cov{r,c} = mat;                   % store thresholded covariance matrix
        clear mask mat
    end
end



 %% Run connectomics on transformed uncoupling matrices
load('roiLabels_wGroupings.mat','roi_labels')
load('colormaps.mat','bluered_cmap')
N = length(roi_labels);
for i=2:rw_transformed
    for j=2:cl_transformed
        matrix = transformed_uncoupled{i,j};
        [agrMat, mrccPartition,allPartitions] = mrcc_wrapper(matrix*-1,10000,1,roi_labels);
        mrccPartition = fcn_relabel_partitions(mrccPartition);
        [X,Y,idx_ord]=grid_communities(mrccPartition);
        f(i,j)= figure('units','inches','position',[1 1 6 6],'paperpositionmode','auto');
        imagesc(matrix(idx_ord,idx_ord)); axis square
        colormap(bluered_cmap); clim([-1 1])
        hold on
        plot(X,Y,'k','LineWidth',2)
        xticks(1:1:N); yticks(1:1:N);
        yticklabels(roi_labels(idx_ord)); xticklabels(roi_labels(idx_ord)); xtickangle(90)
        colorbar
        title(strcat(transformed_uncoupled{i,1}," ",transformed_uncoupled{1,j}," Community Ordered Covariance"))
        exportgraphics(f(i,j),strcat(pwd,"\",transformed_uncoupled{i,1},"_",transformed_uncoupled{1,j},"_coupling.png"),"Resolution",300)
    end
end
close all



save("output.mat","uncoupling_supermatrix","transformed_uncoupled","T1_T2_matrix")
close all
end
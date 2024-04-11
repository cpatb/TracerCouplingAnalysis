function uncoupling_preprocessing(REF,EXP)

[nrA, ncA] = size(REF);
[nrB, ncB] = size(EXP);

%% Perform data checks
% ensure number of regions in arrays matches regions input
% parameter
if nrA ~= nrB
    fprintf(2,'Number of brain regions in experimental group does not match number of brain regions in reference group. Exiting...\n')
    return
end

%% compute the mean, standard deviation of the reference and experimental populations for
% every region
for i=1:nrA
    meanA(i) = mean(REF(i,:)); %#ok<AGROW>
    stdevA(i) = std(REF(i,:)); %#ok<AGROW>
end
for i=1:nrB
    meanB(i) = mean(EXP(i,:)); %#ok<AGROW>
    stdevB(i) = std(EXP(i,:)); %#ok<AGROW>
end
% compute the z-score of within-cohort per-region uptake for ref group, exp group using mean and stdev of ref group 
% ref group
for i=1:nrA
    for j=1:ncA
        z_ref(i,j) = (REF(i,j)-meanA(i))/stdevA(i); %#ok<AGROW>
    end
end

% exp group
for i=1:nrB
    for j=1:ncB
        z_exp2ref(i,j) = (EXP(i,j)-meanA(i))/stdevA(i); %#ok<AGROW>
    end
end

%%
% fetch size of z-score matrices (should be same size as the original
% matrices)
[z_ref_r, z_ref_c] = size(z_ref);
[z_exp2ref_r, z_exp2ref_c] = size(z_exp2ref);

if size(z_ref) ~= size(REF)
    fprintf(2,'Error - Dimensions differ between reference cohort and Z-scored transformation of reference cohort. Exiting...\n')
    return
end

if size(z_exp2ref) ~= size(EXP)
    fprintf(2,'Error - Dimensions differ between experimental cohort and Z-scored transformation of experimental cohort with respect to reference mean and stdev. Exiting...\n')
    return
end


% create vectors of mean, SEM Z-score data
for k=1:z_ref_r
    z_ref_mean(k) = mean(z_ref(k,:)); %#ok<AGROW>
    z_ref_SEM(k) = std(z_ref(k,:))/sqrt(z_ref_c); %#ok<AGROW>
end

for k=1:z_exp2ref_r
    z_exp2ref_mean(k) = mean(z_exp2ref(k,:)); %#ok<AGROW>
    z_exp2ref_SEM(k) = std(z_exp2ref(k,:))/sqrt(z_exp2ref_c); %#ok<AGROW>
end

% compare the z-scores of a specific brain between ref z-scored to ref and exp z-scored to ref
for l=1:z_exp2ref_r
    [h,p,~,~] = ttest2(z_exp2ref(l,:),z_ref(l,:),"Vartype","unequal","Tail","both");
    ttest_of_zscores(l,1) = h; ttest_of_zscores(l,2) = p; %#ok<AGROW>
end
%%
t = array2table(ttest_of_zscores,"RowNames",{'AI','AuDMV','CPu','Cg','CC','DLO','DLIVEnt','DI','ECT','Fornix','FrA','HIP','LO','MO','PtPR','PtA','PRH','PrL','M1','S1','RSC','M2','S2','TeA','TH','VO','V1V2'},"VariableNames",{'Logical Output', 'p-value'});

% compare the FDG and PTSM zscores 

save('preprocessed.mat',"z_ref","z_ref_mean","z_exp2ref","z_exp2ref_mean","z_exp2ref_SEM","t","p","ttest_of_zscores")
end
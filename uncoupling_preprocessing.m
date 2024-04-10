function uncoupling_preprocessing(ref_pop,exp_pop,regions)

[nrA, ncA] = size(ref_pop);
[nrB, ncB] = size(exp_pop);

%% Perform data checks
% ensure number of regions in arrays matches regions input
% parameter
if nrA ~= regions
    fprintf(2,'Number of brain regions in input array A does not match specified region count. Exiting...\n')
    return
end
if nrB ~= regions
    fprintf(2,'Number of brain regions in input array B does not match specified region count. Exiting...\n')
    return
end

%% compute the mean, standard deviation of the reference and experimental populations for
% every region
for i=1:nrA
    meanA(i) = mean(ref_pop(i,:));
    stdevA(i) = std(ref_pop(i,:));
end
for i=1:nrB
    meanB(i) = mean(exp_pop(i,:));
    stdevB(i) = std(exp_pop(i,:));
end
% compute the z-score of within-cohort per-region uptake for ref group, exp group using mean and stdev of ref group 
% ref group
for i=1:nrA
    for j=1:ncA
        z_ref(i,j) = (ref_pop(i,j)-meanA(i))/stdevA(i);
    end
end

% exp group
for i=1:nrB
    for j=1:ncB
        z_exp2ref(i,j) = (exp_pop(i,j)-meanA(i))/stdevA(i);
    end
end

%%
% fetch size of z-score matrices (should be same size as the original
% matrices)
[z_ref_r, z_ref_c] = size(z_ref);
[z_exp2ref_r, z_exp2ref_c] = size(z_exp2ref);

% create vectors of mean, SEM data
for k=1:z_ref_r
    z_ref_mean(k) = mean(z_ref(k,:));
    z_ref_SEM(k) = std(z_ref(k,:))/sqrt(z_ref_c);
end

for k=1:z_exp2ref_r
    z_exp2ref_mean(k) = mean(z_exp2ref(k,:));
    z_exp2ref_SEM(k) = std(z_exp2ref(k,:))/sqrt(z_exp2ref_c);
end

% compare the z-scores of a specific brain between ref z-scored to exp and
% exp z-scored to exp
for l=1:z_exp2ref_r
    [h,p,~,~] = ttest2(z_exp2ref(l,:),z_ref(l,:),"Vartype","unequal","Tail","both");
    ttest_of_zscores(l,1) = h; ttest_of_zscores(l,2) = p;
end
%%
t = array2table(ttest_of_zscores,"RowNames",{'AI','AuDMV','CPu','Cg','CC','DLO','DLIVEnt','DI','ECT','Fornix','FrA','HIP','LO','MO','PtPR','PtA','PRH','PrL','M1','S1','RSC','M2','S2','TeA','TH','VO','V1V2'},"VariableNames",{'Logical Output', 'p-value'});

% compare the FDG and PTSM zscores 

save('preprocessed.mat',"z_ref","z_ref_mean","z_exp2ref","z_exp2ref_mean","z_exp2ref_SEM","t","p","ttest_of_zscores")
end
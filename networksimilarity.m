%% Retrieve two adjacency matrices, sort into positive and negative constituent matrices
covMattransform("covariance_out_tier1_LOAD2 HFD")
load("covMat_constituents.mat")

%% Generate network similarity metrics: Jaccard distance, Euclidean distance
% Jaccard Index 
JaccardPET(ThrPosCovMat,ThrNegCovMat)
load("Jaccard.mat")
thrJaccard = Jaccard; clear Jaccard
JaccardPET(PosCovMat,NegCovMat)
load("Jaccard.mat")

%% Output as Excel File
WithinVarNames = {'Positive Female','Positive Male','Negative Female','Negative Male'};
WithinRowNames =   {'4M/12M','12M/18M','18M/24M','24M/12M HFD','12M HFD/18M HFD'};
BetweenVarNames = {'Positive','Negative'};
BetweenRowNames = {'4M','12M','18M','24M','12M HFD','18M HFD'};
%% Unthresholded
% Within sex sheet generation
w = horzcat(Jaccard{2,2},Jaccard{2,3});
t1 = table(w(:,1),w(:,2),w(:,3),w(:,4),'VariableNames',WithinVarNames,'RowNames',WithinRowNames);
% Between sex sheet generation
b = horzcat(Jaccard{3,2},Jaccard{3,3});
t2 = table(b(:,1),b(:,2),'VariableNames',BetweenVarNames,'RowNames',BetweenRowNames);
% Write out tables
writetable(t1,"Jaccard Indexing.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Unthr Within")
writetable(t2,"Jaccard Indexing.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Unthr Between")
%% Thresholded
% Within sex sheet generation
thrw = horzcat(thrJaccard{2,2},thrJaccard{2,3});
t3 = table(thrw(:,1),thrw(:,2),thrw(:,3),thrw(:,4),'VariableNames',WithinVarNames,'RowNames',WithinRowNames);
% Between sex sheet generation
thrb = horzcat(thrJaccard{3,2},thrJaccard{3,3});
t4 = table(thrb(:,1),thrb(:,2),'VariableNames',BetweenVarNames,'RowNames',BetweenRowNames);
% Write out tables
writetable(t3,"Jaccard Indexing.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Thr Within")
writetable(t4,"Jaccard Indexing.xlsx","FileType","spreadsheet","UseExcel",true,"WriteRowNames",true,"WriteVariableNames",true,"Sheet","Thr Between")
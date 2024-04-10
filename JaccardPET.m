function JaccardPET(PosCovMat,NegCovMat)
%% Jaccard Index 
% Vertically concatenate binarized lower triangle of matrices into a single vector per
% adjacency matrix
[nrP, ncP] = size(PosCovMat);
for i=2:nrP
    for j=2:ncP
        P = PosCovMat{i,j}; catP = P'; catP = catP(:)'; catP = catP.';
        vertcatPos{i,j} = catP;
        N = NegCovMat{i,j}; catN = N'; catN = catN(:)'; catN = catN.';
        vertcatNeg{i,j} = catN;
        for k=1:nrP
            vertcatPos{k,1} = PosCovMat{k,1}; vertcatNeg{k,1} = PosCovMat{k,1};
        end
        for l=1:ncP
            vertcatPos{1,l} = PosCovMat{1,l}; vertcatNeg{1,l} = PosCovMat{1,l};
        end
    end
end
clear nrP ncP P N i j k l
% Compute Jaccard Network Comparison Within, Between Sexes
[nr, nc] = size(vertcatPos);
% Within sex compared to vertically adjacent time points
for i=2:nr-1
    for j=2:nc
        P = vertcatPos{i,j}; N = vertcatNeg{i,j};
        JP(i-1,j-1) = sum(vertcatPos{i,j} & vertcatPos{i+1,j})/sum(vertcatPos{i,j} | vertcatPos{i+1,j});
        JN(i-1,j-1) = sum(vertcatNeg{i,j} & vertcatNeg{i+1,j})/sum(vertcatNeg{i,j} | vertcatNeg{i+1,j});
        Jaccard{2,2} = JP; Jaccard{2,3} = JN;
    end
end
clear i j P N JP JN
% Between sex at a given time
for i=2:nr
    for j=2:nc
        P = vertcatPos{i,j}; N = vertcatNeg{i,j};
        JP(i-1,1) = sum(vertcatPos{i,nc-1} & vertcatPos{i,nc})/sum(vertcatPos{i,nc-1} | vertcatPos{i,nc});
        JN(i-1,1) = sum(vertcatNeg{i,nc-1} & vertcatNeg{i,nc})/sum(vertcatNeg{i,nc-1} | vertcatNeg{i,nc});
        Jaccard{3,2} = JP; Jaccard{3,3} = JN;
    end
end
Jaccard{1,1} = 'Jaccard Index Outputs';
Jaccard{1,2} = 'Positive';
Jaccard{1,3} = 'Negative';
Jaccard{2,1} = 'Within Sex';
Jaccard{3,1} = 'Between Sex';
save("Jaccard.mat","Jaccard")
end
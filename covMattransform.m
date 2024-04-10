function covMattransform(tier1)
%% Load in tier1 unthresholded and thresholded covariance matrix output
load(tier1,"covMat","thr_cov"); A = covMat; T = thr_cov; clear covMat; clear thr_cov;
%% Reformat thresholded covariance matrix output from tier1
[nr,nc] = size(A);
[~,~,nd,ns] = size(T);
for k=1:nd
    for l=1:ns
        a = T(:,:,k,l);
        % line after this is the error - unable to convert from 4D double
        % to cell array
        thr_covMat{k,l} = a; 
    end
end
thr_covMat1=thr_covMat'; thr_covMat2=reshape(thr_covMat1,[nr-1,nc-1]);

[nr,nc] = size(thr_covMat2);
for i=2:nr+1
    for j=2:nc+1
        thr_covMat3{i,j} = thr_covMat2{i-1,j-1};
    end
end
clear nr nc
[nr,nc] = size(A);
for i=1:nr
    thr_covMat3{i,1} = A{i,1};
end
for j=1:nc
    thr_covMat3{1,j} = A{1,j};
end
thr_covMat = thr_covMat3;
%% Transform the unthresholded covariance matrix into two binarized matrices representing positive + negative connections
for i=2:nr
    for j=2:nc
        B = cell2mat(A(i,j));
        ThrB = cell2mat(thr_covMat(i,j));
        [rw,cl] = size(B);
        [rw2,cl2] = size(ThrB);
        if rw==rw2
            if cl==cl2
                for k=1:rw
                    for l=1:cl
                        if B(k,l)>0
                            C(k,l) = 1;
                            C(k,k) = 0;
                        else
                            C(k,l) = 0;
                        end
                        if B(k,l)<0
                            D(k,l) = 1;
                        else
                            D(k,l) = 0;
                        end
                        if ThrB(k,l)>0
                            ThrC(k,l) = 1;
                            ThrC(k,k) = 0;
                        else
                            ThrC(k,l) = 0;
                        end
                        if ThrB(k,l)<0
                            ThrD(k,l) = 1;
                        else
                            ThrD(k,l) = 0;
                        end    
                    end
                end
            else
                fprintf(2,'Number of columns varies. Check that matrix inputs match.\n')
            end
        else
            fprintf(2,'Number of columns varies. Check that matrix inputs match.\n')
        end   
        PosCovMat{i,j} = C; NegCovMat{i,j} = D;
        ThrPosCovMat{i,j} = ThrC; ThrNegCovMat{i,j} = ThrD;
        for m=1:nr
            PosCovMat{m,1} = A{m,1}; NegCovMat{m,1} = A{m,1};
        end
        for n=1:nc
            PosCovMat{1,n} = A{1,n}; NegCovMat{1,n} = A{1,n};
        end
        for m=1:nr
            ThrPosCovMat{m,1} = A{m,1}; ThrNegCovMat{m,1} = A{m,1};
        end
        for n=1:nc
            ThrPosCovMat{1,n} = A{1,n}; ThrNegCovMat{1,n} = A{1,n};
        end
    end
end

%% Save outputs
save('covMat_constituents.mat',"PosCovMat","NegCovMat","ThrPosCovMat","ThrNegCovMat")
end

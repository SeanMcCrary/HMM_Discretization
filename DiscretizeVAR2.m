function [grid,tpm,ll,sigs,Y1,Y2] = DiscretizeVAR2(m,A,Sigma)

%% Arguments 
% m     - number of grid points 
% A     - coefficient of lag matrix 
% Sigma - covariance matrix of errors   

%&------------------------------------------------------------------------%
%% Simulate Data & Construct Initial Guess 
%&------------------------------------------------------------------------%
rng(1992)
N    = 100;
nsim = 1000; % batch t

% ergodic distribution 
sigy2  = reshape((eye(4)-kron(A,A))\Sigma(:),2,2); 
Y0     = mvnrnd([0;0],sigy2,m*nsim);

mu1  = mean(Y0(:,1));
mu2  = mean(Y0(:,2));
std1 = std(Y0(:,1));
std2 = std(Y0(:,2));

D = [(Y0(:,1)-mu1)./std1 (Y0(:,2)-mu2)./std2];
D = max(min(D,3.291),-3.291);

[~, grid] = kmeans(D, m);
grid = [(std1*grid(:,1)+mu1) (std2*grid(:,2)+mu2)];
scatter(Y0(:,1),Y0(:,2)); hold on; scatter(grid(:,1),grid(:,2),'filled'); hold off

%-------------------------------------------------------------------------%
%% Update transition recursively

tol    = 1;
count  = 0;
maxT   = 50000;
grid1  = zeros(m,1);
denom0 = zeros(m,1);
denom1 = zeros(m,1);
num0   = zeros(m,m);
num1   = zeros(m,m);
tpm1   = zeros(m,m);
tpm0   = tpm1;
grid0  = grid;
Y1I    = Y0(1:N,1);
Y10    = [];
Y2I    = Y0(1:N,2);
Y20    = []; 

while tol>10^-2 && count <= maxT

    % Simulate new observations
    Y1      = zeros(N,nsim);
    Y2      = zeros(N,nsim);
    eps    = mvnrnd([0;0],Sigma,N*nsim); 
    epsy1  = reshape(eps(:,1),N,nsim); 
    epsy2  = reshape(eps(:,2),N,nsim); 

    Y1(:,1) = Y1I;
    Y2(:,1) = Y2I;

    for ti=2:nsim
        Y1(:,ti) = A(1,1).*Y1(:,ti-1) + A(1,2).*Y2(:,ti-1) + epsy1(:,ti); 
        Y2(:,ti) = A(2,1).*Y1(:,ti-1) + A(2,2).*Y2(:,ti-1) + epsy2(:,ti); 
    end 

    Y1 = Y1'; Y2=Y2';
    [~,Gidx] = pdist2(grid,[Y1(:) Y2(:)],'euclidean','Smallest',1);
       Gidx  = reshape(Gidx,nsim,N)';
     Y1 = Y1'; Y2=Y2';

    % Compute marginal estimates
    for mi=1:m
        I          = (Gidx==mi);
        denom1(mi) = sum(I,'all');
    end

    for ri=1:m
        for ci=1:m
            I0 = (Gidx(:,1:end-1))==ri;
            if sum(I0) == 0
                num1(ri,ci) = 0;
            else
                num1(ri,ci) = sum(I0.*(Gidx(:,2:end)==ci),'all');
            end
        end
    end
    tpm1 = (num1+num0)./(denom1+denom0);
    tpm1 = tpm1./sum(tpm1,2);

    % Check convergence
    count  = count+nsim;
    tol1   = norm(abs(tpm1-tpm0)./(abs(tpm0)+1));
    tol2   = norm(abs(grid1-grid0)./(abs(grid0)+1));
    tpm0   = tpm1;
    grid0  = grid1;
    denom0 = denom1+denom0;
    num0   = num1+num0;
    Y10    = [Y10 Y1]; %#ok
    Y20    = [Y20 Y2]; %#ok
    Y1I    = Y10(:,end);
    Y2I    = Y20(:,end); 

    printtol = max(tol1,tol2);
    fprintf('Tolerance %4.6f Time T=%4.0f \n',printtol,count)
    tol      = printtol;

end

% Initial for EM 
[N,T]  = size(Y10);

sigY1   = sqrt(mean((Y1(:)-grid(Gidx,1)).^2,'all')).*ones(m,T);
sigY2   = sqrt(mean((Y2(:)-grid(Gidx,2)).^2,'all')).*ones(m,T);

[~,Idx] = sort(grid(:,1)); 
gridy1  = repmat(grid(Idx,1),1,T);
gridy2  = repmat(grid(Idx,2),1,T);
Pi      = tpm1(Idx,Idx); 

pistat = ones(1,m)/(eye(m)-Pi+ones(m,m)); 
Y1     = Y10';
Y2     = Y20'; 

tpm = cell(T-1,1);
for ti = 1:T-1
    tpm{ti} = Pi;
end

%&------------------------------------------------------------------------%
%%  EM Algorithm 
%&------------------------------------------------------------------------%

% Forward and backward probs
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
cscl    = cell(T,1);

maxiter       = 1000;
loglikelihood = zeros(maxiter,1);

tol     = 1;

for niter = 1:maxiter

    tic
    %% forward
    sigs      = sigY1(:,1);
    phiy      = (exp(-((Y1(1,:)-gridy1(:,1)).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    sigs      = sigY2(:,1);
    phiz      = (exp(-((Y2(1,:)-gridy2(:,1)).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    phi{1}    = phiy.*phiz;

    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs       = sigY1(:,ti);
        phiy       = (exp(-((Y1(ti,:)-gridy1(:,ti)).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        sigs       = sigY2(:,ti);
        phiz       = (exp(-((Y2(ti,:)-gridy2(:,ti)).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        phi{ti}    = phiy.*phiz;

        pis        = tpm{ti-1};
        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},1e-7);
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = alf{ti}.*cscl{ti};

        llk        = -log(cscl{ti});
        ll         = ll + sum(llk,"all");

    end

    %% backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};

    gtemp   = alfscl{T}.*betscl{T};
    gam{T}  = gtemp./sum(gtemp,1);

    for ti = (T-1):-1:1

        pis        = tpm{ti};
        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),1e-7);
        betscl{ti} = bet{ti}.*cscl{ti};

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

    end
    ittoc = toc;
   
    %% Update Mean (i.e. Grid)
    munum   = zeros(m,1); 
    mudenom = zeros(m,1); 
    for ti=1:T
        munum      = munum + sum(gam{ti}.*Y1(ti,:),2); 
        mudenom    = mudenom + sum(gam{ti},2);
    end 
    grid0  = munum./mudenom;
    gridy1 = repmat(grid0,1,T);

    munum   = zeros(m,1);
    mudenom = zeros(m,1);
    for ti=1:T
        munum      = munum + sum(gam{ti}.*Y2(ti,:),2);
        mudenom    = mudenom + sum(gam{ti},2);
    end
    grid0  = munum./mudenom;
    gridy2 = repmat(grid0,1,T);

    %% Update Transitions
    tpmnum  = zeros(m,m);
    for ti=1:(T-1)
        xis = alfscl{ti} * (phi{ti+1} .* betscl{ti+1})' .* Pi;
        tpmnum  = tpmnum + xis;
    end
    Pi     = tpmnum./sum(tpmnum,2);
    for ti=1:T-1
        tpm{ti} = Pi;
    end
    pistat = ones(1,m)/(eye(m)-Pi+ones(m,m));  % update initial

    %% Update sigY and sigZ  
    sse        = cell(T,1);
    ssecumul   = zeros(m,N); 
    wgtcumul   = zeros(m,N); 
    for ti = 1:T
        sse{ti}  = gam{ti}.*((Y1(ti,:)-gridy1(:,ti)).^2);
        ssecumul = ssecumul + sse{ti}; 
        wgtcumul = wgtcumul + gam{ti}; 
    end 
    sigY1 = sqrt(sum(ssecumul,'all')./sum(wgtcumul,'all')).*ones(m,T); 

    sse        = cell(T,1);
    ssecumul   = zeros(m,N); 
    wgtcumul   = zeros(m,N); 
    for ti = 1:T
        sse{ti}  = gam{ti}.*((Y2(ti,:)-gridy2(:,ti)).^2);
        ssecumul = ssecumul + sse{ti}; 
        wgtcumul = wgtcumul + gam{ti}; 
    end 
    sigY2 = sqrt(sum(ssecumul,'all')./sum(wgtcumul,'all')).*ones(m,T);

    % Print
    loglikelihood(niter,1) = ll;

    if niter > 2
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end

    fprintf('Iteration %4.0f Update Error %4.6f  Time %4.2f LogLik %4.6f \n',niter,tol,ittoc,ll)

    if tol<10^-3
        break
    end

end

%% output 
tpm    = Pi; 
grid   = [gridy1(:,1) gridy2(:,1)]; 
sigs   = [sigY1(1,1)  sigY2(1,1)]; 

end 

function [out] = DiscretizeGKOS(m1,m2)

% m1: number of 0 earnings states 
% m2: number of postive earnings states 

% check dimension
m = m1 + m2;
if m>50
    fprintf('Warning! This many grid points will require more simulated data! \n')
    out = []; 
    return  
end

fprintf('Initializing... \n')

load GKOS_SIM %#ok 

y   = log(y'+1);
z   = z';
v   = v';
yu  = log(exp(z)+1);  % y without unemployment shock

% Some trimming
[N, T] = size(y);

eta   = zeros(N,T);
for ti=1:T
    lbeta     = 0;
    ubeta     = mean(y(:,ti)) + norminv(0.9995)*std(y(:,ti));
    eta(:,ti) = max(min(y(:,ti),ubeta),lbeta);
end

y    = eta;

%-------------------------------------------------------------------------%
%  Parameters and initial guess
maxiter = 1000;

% Quantile binning
pistat  = ones(1,m)./m;

% bin the non-zero values of y by quantile
q       = linspace(0,1,m2+1);
bin     = zeros(m2+1,T);
for ti=1:T
    bin(:,ti) = quantile(y(y(:,ti)>0,ti),q);
end
bin(1,:) = 1e-8;

numzero = m1*ones(1,T);

grid  = zeros(m,T);
gridz = zeros(m,T);

for ti = 1:T
    for mi=1:m2
        I              = (y(:,ti)>=bin(mi,ti)).*(y(:,ti)<=bin(mi+1,ti));
        grid(m1+mi,ti) =  mean(y(I>0,ti));
    end
end

bin0 = cell(1,T);
for ti=1:T
    bin0{ti} = quantile(yu(v(:,ti)==1,ti),linspace(0,1,numzero(ti)+1));
end

for ti = 1:T

    for mi=1:m1
        I            = (v(:,ti)==1).*(yu(:,ti)>=bin0{ti}(mi)).*(yu(:,ti)<=bin0{ti}(mi+1));
        gridz(mi,ti) =  mean(z(I>0,ti));
        if ti==1
            pistat(1,mi) = sum(I)/N;
        end
    end

    for mi=1:m2
        I               = (y(:,ti)>=bin(mi,ti)).*(y(:,ti)<=bin(mi+1,ti));
        gridz(m1+mi,ti) =  mean(z(I>0,ti));
        if ti==1
            pistat(1,m1+mi) = sum(I)/N;
        end
    end

end

% transition matrix
tpm  = cell(1,T-1);
for ti = 1:T-1
    tpm{ti} = zeros(m,m);
end

for ti = 1:T-1
    for ri=1:m
        for ci=1:m

            if grid(ri,ti)>0
                I0 =  (y(:,ti)>bin(ri-m1,ti)).*(y(:,ti)<=bin(ri+1-m1,ti));
            else
                I0 = (v(:,ti)==1).*(yu(:,ti)>=bin0{ti}(ri)).*(yu(:,ti)<=bin0{ti}(ri+1));
            end

            if sum(I0) == 0
                tpm{ti}(ri,ci) = 0;
            else

                if grid(ci,ti+1)>0
                    I1 =  (y(:,ti+1)>bin(ci-m1,ti+1)).*(y(:,ti+1)<=bin(ci-m1+1,ti+1));
                else
                    I1 = (v(:,ti+1)==1).*(yu(:,ti+1)>=bin0{ti+1}(ci)).*(yu(:,ti+1)<=bin0{ti+1}(ci+1));
                end

                tpm{ti}(ri,ci) = sum(I1.*I0)/sum(I0);

            end
        end
    end

    tpm{ti} = tpm{ti}./sum(tpm{ti},2);

end

Y      = eta';
Z      = z';

sigY   = (std(Y(:))/m).*ones(m,T);
sigZ   = (std(Z(:))/m).*ones(m,T);

%-------------------------------------------------------------------------%
%% EM Algorithm - Binning
%-------------------------------------------------------------------------%

% Forward and backward probs
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
cscl    = cell(T,1);

loglikelihood=zeros(maxiter,1);

tol =1;

for niter = 1:maxiter

    %% forward
    sigs      = sigY(:,1);
    gs        = grid(:,1);
    phiy      = (exp(-((Y(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);

    sigs      = sigZ(:,1);
    gs        = gridz(:,1);
    phiz      = (exp(-((Z(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    phi{1}    = phiy.*phiz;

    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs      = sigY(:,ti);
        gs        = grid(:,ti);
        phiy      = (exp(-((Y(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);

        sigs      = sigZ(:,ti);
        gs        = gridz(:,ti);
        phiz      = (exp(-((Z(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        phi{ti}   = phiy.*phiz;

        pis        = tpm{ti-1};
        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},10^-7);
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = alf{ti}.*cscl{ti};

        llk       = -log(cscl{ti});
        ll        = sum(llk,"all");

    end

    %% backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};

    gtemp      = alfscl{T}.*betscl{T};
    gam{T}     = gtemp./sum(gtemp,1);

    gs               = grid(:,T);
    sigY(1:m1,T)     = sqrt(sum(gam{T}(1:m1,:).*((Y(T,:)-gs(1:m1)).^2),'all')./(sum(gam{T}(1:m1,:),'all')));
    sigY((m1+1):m,T) = sqrt(sum(gam{T}((m1+1):m,:).*((Y(T,:)-gs((m1+1):m)).^2),'all')./(sum(gam{T}((m1+1):m,:),'all')));

    gs               = gridz(:,T);
    sigZ(1:m1,T)     = sqrt(sum(gam{T}(1:m1,:).*((Z(T,:)-gs(1:m1)).^2),'all')./(sum(gam{T}(1:m1,:),'all')));
    sigZ((m1+1):m,T) = sqrt(sum(gam{T}((m1+1):m,:).*((Z(T,:)-gs((m1+1):m)).^2),'all')./(sum(gam{T}((m1+1):m,:),'all')));

    for ti = (T-1):-1:1

        pis        = tpm{ti};
        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),10^-7);
        betscl{ti} = bet{ti}.*cscl{ti};

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

        gs                = grid(:,ti);
        sigY(1:m1,ti)     = sqrt(sum(gam{ti}(1:m1,:).*((Y(ti,:)-gs(1:m1)).^2),'all')./(sum(gam{ti}(1:m1,:),'all')));
        sigY((m1+1):m,ti) = sqrt(sum(gam{ti}((m1+1):m,:).*((Y(ti,:)-gs((m1+1):m)).^2),'all')./(sum(gam{ti}((m1+1):m,:),'all')));

        gs                = gridz(:,ti);
        sigZ(1:m1,ti)     = sqrt(sum(gam{ti}(1:m1,:).*((Z(ti,:)-gs(1:m1)).^2),'all')./(sum(gam{ti}(1:m1,:),'all')));
        sigZ((m1+1):m,ti) = sqrt(sum(gam{ti}((m1+1):m,:).*((Z(ti,:)-gs((m1+1):m)).^2),'all')./(sum(gam{ti}((m1+1):m,:),'all')));

    end

    % Print
    loglikelihood(niter,1) = ll;

    if niter > 9

        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);

    end

    if mod(niter,5)==0
        fprintf('Bin-EM Iteration %4.0f Update Error %4.6f LogLik %4.6f \n',niter,tol,ll)
    end 

    if tol<10^-4
        break
    end

end

%-------------------------------------------------------------------------%
% EM Algorithm - HMM
%-------------------------------------------------------------------------%

% Forward and backward probs
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
muyt    = zeros(m,T);
muzt    = zeros(m,T);
cscl    = cell(T,1);

loglikelihood=zeros(maxiter,1);

tol =1;

for niter = 1:maxiter

    % forward
    sigs      = sigY(:,1);
    gs        = grid(:,1);
    phiy      = (exp(-((Y(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);

    sigs      = sigZ(:,1);
    gs        = gridz(:,1);
    phiz      = (exp(-((Z(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);

    phi{1}    = phiy.*phiz;

    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs      = sigY(:,ti);
        gs        = grid(:,ti);
        phiy      = (exp(-((Y(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);

        sigs      = sigZ(:,ti);
        gs        = gridz(:,ti);
        phiz      = (exp(-((Z(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        phi{ti}   = phiy.*phiz;

        pis        = tpm{ti-1};
        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},10^-7);
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = alf{ti}.*cscl{ti};

        llk       = -log(cscl{ti});
        ll        = sum(llk,"all");

    end

    % backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};

    gtemp      = alfscl{T}.*betscl{T};
    gam{T}     = gtemp./sum(gtemp,1);

    muyt(:,T)        = sum(gam{T}.*Y(T,:),2)./sum(gam{T},2);
    gs               = grid(:,T);
    sigY(1:m1,T)     = sqrt(sum(gam{T}(1:m1,:).*((Y(T,:)-gs(1:m1)).^2),'all')./(sum(gam{T}(1:m1,:),'all')));
    sigY((m1+1):m,T) = sqrt(sum(gam{T}((m1+1):m,:).*((Y(T,:)-gs((m1+1):m)).^2),'all')./(sum(gam{T}((m1+1):m,:),'all')));

    muzt(:,T)        = sum(gam{T}.*Z(T,:),2)./sum(gam{T},2);
    gs               = gridz(:,T);
    sigZ(1:m1,T)     = sqrt(sum(gam{T}(1:m1,:).*((Z(T,:)-gs(1:m1)).^2),'all')./(sum(gam{T}(1:m1,:),'all')));
    sigZ((m1+1):m,T) = sqrt(sum(gam{T}((m1+1):m,:).*((Z(T,:)-gs((m1+1):m)).^2),'all')./(sum(gam{T}((m1+1):m,:),'all')));

    for ti = (T-1):-1:1

        pis        = tpm{ti};
        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),10^-7);
        betscl{ti} = bet{ti}.*cscl{ti};

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

        muyt(:,ti)        = sum(gam{ti}.*Y(ti,:),2)./sum(gam{ti},2);
        gs                = grid(:,ti);
        sigY(1:m1,ti)     = sqrt(sum(gam{ti}(1:m1,:).*((Y(ti,:)-gs(1:m1)).^2),'all')./(sum(gam{ti}(1:m1,:),'all')));
        sigY((m1+1):m,ti) = sqrt(sum(gam{ti}((m1+1):m,:).*((Y(ti,:)-gs((m1+1):m)).^2),'all')./(sum(gam{ti}((m1+1):m,:),'all')));

        muzt(:,ti)        = sum(gam{ti}.*Z(ti,:),2)./sum(gam{ti},2);
        gs                = gridz(:,ti);
        sigZ(1:m1,ti)     = sqrt(sum(gam{ti}(1:m1,:).*((Z(ti,:)-gs(1:m1)).^2),'all')./(sum(gam{ti}(1:m1,:),'all')));
        sigZ((m1+1):m,ti) = sqrt(sum(gam{ti}((m1+1):m,:).*((Z(ti,:)-gs((m1+1):m)).^2),'all')./(sum(gam{ti}((m1+1):m,:),'all')));

        % transitions
        tpmnum  = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*tpm{ti};
        tpm{ti} = tpmnum./sum(tpmnum,2);

    end

    % stationary distribution
    pistat = mean(gam{1},2)';

    % Update Mean (i.e. Grid)
    for ti=1:T
        muyt(:,ti)        = sort(muyt(:,ti));
        muzt((m1+1):m,ti) = sort(muzt((m1+1):m,ti));
        muzt(1:m1,ti)     = sort(muzt(1:m1,ti));
    end
    for mi=1:m
        muyt(mi,:) = smooth(muyt(mi,:),3);
        muzt(mi,:) = smooth(muzt(mi,:),3);
    end
    grid(m1+1:end,:)  = muyt(m1+1:end,:);
    gridz = muzt;

    % Print
    if mod(niter,5) ==0 
        fprintf('HMM-EM Iteration %4.0f Update Error %4.6f LogLik %4.6f \n',niter,tol,ll)
    end 
    
    loglikelihood(niter,1) = ll;

    if niter > 9
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end

    if tol<10^-3
        break
    end

end

fprintf('HMM Converged \n')

% output 
out.Grid    = grid; 
out.Pi      = tpm; 
out.initial = pistat; 

end
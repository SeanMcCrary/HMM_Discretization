function [out] = DiscretizeABB(m)

load ABB_SIM %#ok 

% Some trimming
[N, T] = size(eta0); %#ok 

mueta  = mean(eta0);
stdeta = std(eta0);

eta   = zeros(N,T);
for ti=1:T
    lbeta     = mueta(ti) - norminv(0.9995).*stdeta(ti);
    ubeta     = mueta(ti) + norminv(0.9995).*stdeta(ti);
    eta(:,ti) = max(min(eta0(:,ti),ubeta),lbeta);
end
stdeta = std(eta);

fprintf('Initializing... \n')

%&------------------------------------------------------------------------%
%  Parameters and initial guess
%&------------------------------------------------------------------------%

% binning
grid = zeros(m,T);
bins = zeros(m+1,T);

for ti=1:T
    bins(:,ti) = linspace(quantile(eta(:,ti),0.0001),quantile(eta(:,ti),0.9999),m+1);
end

pistat = zeros(1,m);
for ti = 1:T
    for mi=1:m
        I           = (eta(:,ti)>=bins(mi,ti)).*(eta(:,ti)<=bins(mi+1,ti));
        grid(mi,ti) =  mean(eta(I>0,ti));
        if ti==1
            pistat(mi)=mean(I);
        end
    end
end

tpm = cell(1,T-1);

for ti = 1:T-1
    tpm{ti} = zeros(m,m);
end

for ti = 1:T-1
    for ri=1:m
        for ci=1:m
            I0 = (eta(:,ti)>=bins(ri,ti)).*(eta(:,ti)<=bins(ri+1,ti));
            if sum(I0) == 0
                tpm{ti}(ri,ci) = 0;
            else
                I1 = (eta(:,ti+1)>=bins(ci,ti+1)).*(eta(:,ti+1)<=bins(ci+1,ti+1));
                tpm{ti}(ri,ci) = sum(I1.*I0)/sum(I0);
            end
        end
    end
    tpm{ti} = tpm{ti}./sum(tpm{ti},2);
end

Y      = eta';

sigV   = zeros(m,T);
for ti = 1:T
    sigV(:,ti) = m*stdeta(ti);
end

%&------------------------------------------------------------------------%
%   EM Algorithm - Binning
%&------------------------------------------------------------------------%

maxiter = 1000; 

% Forward and backward probs
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
mut     = zeros(m,T);
cscl    = cell(T,1);

loglikelihood = zeros(maxiter,1);

tol     = 1;

for niter = 1:maxiter

    %% forward
    sigs      = sigV(:,1);
    gs        = grid(:,1);
    phi{1}    = (exp(-((Y(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs       = sigV(:,ti);
        gs         = grid(:,ti);
        pis        = tpm{ti-1};
        phi{ti}    = (exp(-((Y(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},10^-7);
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = alf{ti}.*cscl{ti};

        llk       = -log(cscl{ti});
        ll        = sum(llk,"all");

    end

    %% backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};

    gtemp   = alfscl{T}.*betscl{T};
    gam{T}  = gtemp./sum(gtemp,1);

    mut(:,T)   = sum(gam{T}.*Y(T,:),2)./sum(gam{T},2);
    gs         = mut(:,T);
    sigV(:,T)  = sqrt(sum(gam{T}.*((Y(T,:)-gs).^2),'all')./(sum(gam{T},'all')));

    for ti = (T-1):-1:1

        pis        = tpm{ti};
        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),10^-7);
        betscl{ti} = bet{ti}.*cscl{ti};

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

        gs         = grid(:,ti);
        sigV(:,ti) = sqrt(sum(gam{ti}.*((Y(ti,:)-gs).^2),'all')./(sum(gam{ti},'all')));

    end

    pistat = mean(gam{1},2)';
    pistat = pistat./sum(pistat);

    % Print
    loglikelihood(niter,1) = ll;

    if niter > 2
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end

    if mod(niter,5)==0
        fprintf('Bin-EM Iteration %4.0f Update Error %4.6f LogLik %4.6f \n',niter,tol,ll)
    end 

    if tol<10^-3
        break
    end

end

%&------------------------------------------------------------------------%
%   EM Algorithm - JM
%&------------------------------------------------------------------------%

% Forward and backward probs
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
mut     = zeros(m,T);
cscl    = cell(T,1);

loglikelihood = zeros(maxiter,1);

tol     = 1;

for niter = 1:maxiter

    %% forward
    sigs      = sigV(:,1);
    gs        = grid(:,1);
    phi{1}    = (exp(-((Y(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs       = sigV(:,ti);
        gs         = grid(:,ti);
        pis        = tpm{ti-1};
        phi{ti}    = (exp(-((Y(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        alf{ti}    = ((pis')*alfscl{ti-1}).*phi{ti};
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = max(alf{ti}.*cscl{ti},1e-7);

        llk       = -log(cscl{ti});
        ll        = sum(llk,"all");

    end

    %% backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};

    gtemp   = alfscl{T}.*betscl{T};
    gam{T}  = gtemp./sum(gtemp,1);

    mut(:,T)   = sum(gam{T}.*Y(T,:),2)./sum(gam{T},2);
    gs         = mut(:,T);
    sigV(:,T)  = sqrt(sum(gam{T}.*((Y(T,:)-gs).^2),'all')./(sum(gam{T},'all')));

    for ti = (T-1):-1:1

        pis        = tpm{ti};
        bet{ti}    = pis*(betscl{ti+1}.*phi{ti+1});
        betscl{ti} = max(bet{ti}.*cscl{ti},1e-7);

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

        mut(:,ti)  = sort(sum(gam{ti}.*Y(ti,:),2)./sum(gam{ti},2));
        gs         = grid(:,ti);
        sigV(:,ti) = sqrt(sum(gam{ti}.*((Y(ti,:)-gs).^2),'all')./(sum(gam{ti},'all')));

        % transitions
        tpmnum  = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*tpm{ti};
        tpm{ti} = tpmnum./sum(tpmnum,2);

    end

    % update initial distribution
    pistat = mean(gam{1},2)';
    pistat = pistat./sum(pistat);

    % Update Mean (i.e. Grid)
    mut2  = zeros(m,T);
    for ti=1:T
        mut2(:,ti) = sort(mut(:,ti));
    end
    for mi=1:m
        mut2(mi,:) = smooth(mut2(mi,:),3);
    end
    grid = mut2;


    % Print
    loglikelihood(niter,1) = ll;

    if niter > 2
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
        if loglikelihood(niter,1)<loglikelihood(niter-1,1)
            break
        end
    end

    % Print
    if mod(niter,5) ==0 
        fprintf('HMM-EM Iteration %4.0f Update Error %4.6f LogLik %4.6f \n',niter,tol,ll)
    end 

    if tol< 10^-3
        break
    end

end

fprintf('HMM Converged \n')

% output
out.Grid    = grid;
out.Pi      = tpm;
out.initial = pistat;

end
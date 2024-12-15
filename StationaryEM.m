function [Pi,Grid,ll,sige] = StationaryEM(Ysim,Pi0,Grid0)

%% Arguments
% Ysim:  an T-N-k array of data
% Pi0:   an m-m initial guess of transitions
% Grid0: an m-k initial guess of grid points

%% Output
% Pi:   transition probability matrix
% Grid: matrix of estimated grid points
% ll:   log-likelihood
% sige: vector of standard deviations of the error terms

%_________________________________________________________________________%
%% Preliminaries
[T,N,k]  = size(Ysim);
[m,km]   = size(Grid0);
[mr,mc]  = size(Pi0);

if km~=k
    frpritf('ERROR: Dimensions of grid and data do not agree')
    Pi=[]; Grid=[]; ll=[]; sige=[];
    return
elseif m~=mr
    frpritf('ERROR: Dimensions of grid and transition matrix do not agree')
    Pi=[]; Grid=[]; ll=[]; sige=[];
    return
elseif mr~=mc
    frpritf('ERROR: Transition matrix must be square')
    Pi=[]; Grid=[]; ll=[]; sige=[];
    return
elseif T*N < 10000
    fprintf('WARNING: You may want to increase the simulated sample size for more reliable results')
end

% simulated data
Y = cell(k,1);
for ki=1:k
    Y{ki} = squeeze(Ysim(:,:,ki));
end

% initialize sigma
D = [];
for ki=1:k
    Ytemp = Y{ki}';
    D = [D Ytemp(:)]; %#ok
end
[~,Gidx] = pdist2(Grid0,D,'euclidean','Smallest',1);
Gidx  = reshape(Gidx,T,N)';
sige = zeros(m,k);
for ki=1:k
    sige(:,ki)=sqrt(mean((D(:,ki)-Grid0(Gidx,ki)).^2,'all'));
end

% remaining initials
pis    = Pi0;
grid   = Grid0;
pistat = ones(1,m)/(eye(m)-pis+ones(m,m));

%%_________________________________________________________________________%
%&
%&  EM Algorithm
%&_________________________________________________________________________%

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
tol           = 1;

for niter = 1:maxiter
    %__________________________________________________________________________%
    %   E-Step

    tic
    % forward
    phi{1}   = ones(m,N);
    for ki=1:k
        sigs   = sige(:,ki);
        D      = Y{ki};
        gs     = grid(:,ki);
        phiD   = (exp(-((D(1,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        phi{1} = phi{1}.*phiD;
    end

    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        phi{ti}   = ones(m,N);
        for ki=1:k
            sigs    = sige(:,ki);
            D       = Y{ki};
            gs      = grid(:,ki);
            phiD    = (exp(-((D(ti,:)-gs).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
            phi{ti} = phi{ti}.*phiD;
        end

        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},10^-7);
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = alf{ti}.*cscl{ti};
        llk        = -log(cscl{ti});
        ll         = ll + sum(llk,"all");

    end

    % backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};
    gtemp      = alfscl{T}.*betscl{T};
    gam{T}     = gtemp./sum(gtemp,1);

    for ti = (T-1):-1:1

        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),10^-7);
        betscl{ti} = bet{ti}.*cscl{ti};
        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

    end
    ittoc = toc;

    %______________________________________________________________________%
    %   M-Step

    % Update Mean (i.e. Grid)
    munum   = zeros(m,k);
    mudenom = zeros(m,1);
    for ti=1:T
        for ki=1:k
            munum(:,ki) = munum(:,ki) + sum(gam{ti}.*Y{ki}(ti,:),2);
        end
        mudenom    = mudenom + sum(gam{ti},2);
    end
    grid  = munum./mudenom;

    % Update transition probability matrix
    tpmnum  = zeros(m,m);
    for ti=1:(T-1)
        xis = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*pis;
        tpmnum  = tpmnum + xis;
    end
    pis    = tpmnum./sum(tpmnum,2);
    pistat = ones(1,m)/(eye(m)-pis+ones(m,m)); % update initial

    % Update sigma_e
    sse = zeros(m,k);
    for ti = 1:T
        for ki=1:k
            ssetemp   = sum(gam{ti}.*((Y{ki}(ti,:)-grid(:,ki)).^2),'all');
            sse(:,ki) = sse(:,ki) + ssetemp;
        end
    end
    sige = sqrt(sse./(N*T));

    % Check for convergence
    loglikelihood(niter,1) = ll;
    if niter > 1
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end
    fprintf('Iteration %4.0f Update Error %4.6f  Time %4.2f LogLik %4.6f \n',niter,tol,ittoc,ll)
    if tol<10^-4
        break
    end

end

% Output
ll   = ll/(N*T);
Pi   = pis;
Grid = grid;

end
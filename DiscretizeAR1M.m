function [grid,tpm,ll,sigs,Y] = DiscretizeAR1M(m,rho,mixp,mixmu,mixsig)

%% Arguments 
% m      - number of grid points 
% rho    - persistence 
% mixp   - vector of mixture weights 
% mixmu  - vector of mixture means 
% mixsig - vector of mixture standard deviations 

%% Output 
% grid - vector of means 
% tpm  - transition probability matrix
% ll   - loglieklihood 
% sigs - standard deviation of the error term 
% Y    - simulated data set 

%&------------------------------------------------------------------------%
%% Simulate Data & Initial Guess 
%&------------------------------------------------------------------------%

% Grid 
k      = length(mixp);
sige2  = sum(mixp.*(mixsig.^2 + mixmu.^2)) - sum(mixp.*mixmu).^2;
sigy   = sqrt(sige2/(1-rho^2)); 
space  = min(4,1.2*log(m-1));
grid   = linspace(-space*sigy,space*sigy,m)';

% Transition 
Pi = zeros(m,m);
mu = grid; 
for ci = 1:m
    for ri=1:m
        p          = 0;
        for ki=1:k
            p = p + mixp(ki)*normpdf(mu(ci),rho*mu(ri)+mixmu(ki),mixsig(ki));
        end 
        Pi(ri,ci) = p; 
    end 
end  
Pi = Pi./sum(Pi,2); 

%_______________________________________________________________________%
% Initial value of Y 
rng(1992)
N       = 100; 
burn    = 1000; 
Y       = zeros(N,1+burn); 
Y(:,1)  = normrnd(0,sigy,[N 1]); 

eps         = cell(k,1); 
for ki=1:k
    eps{ki} = normrnd(mixmu(ki),mixsig(ki),[N (1+burn)]);
end

for ti=2:(1+burn)
    X       = mnrnd(1,mixp,N); 
    e       = zeros(N,1); 
    for ki=1:k 
        e  = e + X(:,ki).*eps{ki}(:,ti);
    end 
    Y(:,ti) = rho.*Y(:,ti-1) + e; 
end 

YI = Y(:,end); 

%_________________________________________________________________________%
% Recursive binning

nsim = 1000; % batch t
tol   = 1;
count = 0;
maxT  = 30000;

% pre-allocate 
denom0 = zeros(m,1);
denom1 = zeros(m,1);
num0   = zeros(m,m);
num1   = zeros(m,m);
tpm0   = zeros(m,m);
grid0  = grid;
grid1  = grid; 
Y0     = []; 

%_________________________________________________________________________%
% Update transitions recursively

while tol>10^-2 && count <= maxT

    % Simulate new observations
    Y      = zeros(N,nsim+1);
    Y(:,1) = YI; 

    eps         = cell(k,1);
    for ki=1:k
        eps{ki} = normrnd(mixmu(ki),mixsig(ki),[N (1+nsim)]);
    end

    for ti=2:(1+nsim)
        X       = mnrnd(1,mixp,N);
        e       = zeros(N,1);
        for ki=1:k
            e  = e + X(:,ki).*eps{ki}(:,ti);
        end
        Y(:,ti) = rho.*Y(:,ti-1) + e;
    end
    Y        = Y(:,2:end); 
    Y        = Y';
    [~,Yidx] = pdist2(grid,Y(:),'euclidean','Smallest',1);
    Yidx     = reshape(Yidx,nsim,N)';
    Y        = Y';

    % Compute marginal estimates
    for ri=1:m
        I0         = (Yidx(:,1:end-1)==ri);
        denom1(ri) = sum(I0,'all');
        for ci=1:m
            if denom1(ri) == 0
                num1(ri,ci) = 0;
            else
                num1(ri,ci) = sum(I0.*(Yidx(:,2:end)==ci),'all');
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
    Y0     = [Y0 Y]; %#ok
    YI     = Y0(:,end);

    printtol = max(tol1,tol2);
    fprintf('Tolerance %4.6f Time T=%4.0f \n',printtol,count)
    tol      = printtol; 

end

[N, T]    = size(Y0); 
grid      = grid(:);
pis       = Pi; 
pistat    = max(ones(1,m)/(eye(m)-Pi+ones(m,m)),1e-7);
sigs      = ones(m,1)*sqrt(mean((Y-grid(Yidx)).^2,'all'));
Y         = Y0';

tpm = cell(1,T-1);
for ti = 1:T-1
    tpm{ti} = Pi;
end

%&_________________________________________________________________________%
%%
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

    tic
    %% forward
    phi{1}    = (exp(-((Y(1,:)-grid).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        phi{ti}    = (exp(-((Y(ti,:)-grid).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        alf{ti}    = ((pis')*alfscl{ti-1}).*phi{ti};
        cscl{ti}   = 1./sum(alf{ti},1);
        alfscl{ti} = max(alf{ti}.*cscl{ti},1e-7);
        llk        = -log(cscl{ti});
        ll         = ll + sum(llk,"all");

    end

    %% backward
    bet{T}     = ones(m,N);
    betscl{T}  = bet{T}.*cscl{T};
    gtemp      = alfscl{T}.*betscl{T};
    gam{T}     = gtemp./sum(gtemp,1);

    for ti = (T-1):-1:1

        bet{ti}    = pis*(betscl{ti+1}.*phi{ti+1});
        betscl{ti} = max(bet{ti}.*cscl{ti},1e-7);
        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

    end
    ittoc = toc;
   
    % Update mean (i.e. grid)
    munum   = zeros(m,1); 
    mudenom = zeros(m,1); 
    for ti=1:T
        munum      = munum + sum(gam{ti}.*Y(ti,:),2); 
        mudenom    = mudenom + sum(gam{ti},2);
    end 
    grid = sort(munum./mudenom);

    % Update transition probability matrix 
    tpmnum  = zeros(m,m);
    for ti=1:(T-1)
        xis = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*pis;
        tpmnum  = tpmnum + xis;
    end
    pis    = tpmnum./sum(tpmnum,2);
    pistat = max(ones(1,m)/(eye(m)-Pi+ones(m,m)),1e-7); % update initial

    % Update sigma
    ssecumul = 0;
    for ti = 1:T
        ssetemp  = sum(gam{ti}.*((Y(ti,:)-grid).^2),'all');
        ssecumul = ssecumul + ssetemp;
    end
    sigs = sqrt(sum(ssecumul,'all')./(N*T)).*ones(m,1);

    % Check convergence 
    loglikelihood(niter,1) = ll;
    if niter > 2
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end
    fprintf('Iteration %4.0f Update Error %4.6f  Time %4.2f LogLik %4.6f \n',niter,tol,ittoc,ll)
    if tol<0.005
        break
    end

end

% output 
grid = grid(:); 
tpm  = pis; 
sigs = sigs(1); 
ll   = ll/(N*T); 

end 

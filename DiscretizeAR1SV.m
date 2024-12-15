function [gridz,gridy,tpm,ll,sigs,Y,Z] = DiscretizeAR1SV(m,rhoM,rhoV,omega,eta)

%% Arguments 
% m     - number of grid points 
% rhoM  - persistence of levels 
% rhoV  - persistence of volatility
% omega - s.d. of volatiltiy shocks
% eta   - mean volatility
% N     - number of series 
% T     - length of series 

%% Simulation loop
N    = 100;
nsim = 1000; % batch t

%-------------------------------------------------------------------%
%  K-Means on ergodic set 

% Simulate observations
rng(1992)
T      = m*nsim;
burnin = ceil(100./(1-max(rhoM,rhoV)));
erg    = 100; % space between dependent observation for ergodic set
Y      = zeros(1,burnin+erg*T);
Z      = zeros(1,burnin+erg*T);
epsy   = normrnd(0,1,[1 (burnin+erg*T)]);
epsz   = normrnd(0,1,[1 (burnin+erg*T)]);

Z(:,1) = (1-rhoV)*eta + omega.*epsz(:,1);
Y(:,1) = sqrt(Z(:,1)).*epsy(:,1);

for ti=2:(burnin+erg*T)
    Z(:,ti) = (1-rhoV)*eta + rhoV.*Z(:,ti-1) + omega.*epsz(:,ti);
    Y(:,ti) = rhoM.*Y(:,ti-1) + sqrt(Z(:,ti)).*epsy(:,ti);
end

Y    = Y(:,(burnin+1):erg:end);
Z    = Z(:,(burnin+1):erg:end);
muy  = mean(Y(:));
muz  = mean(Z(:));
stdy = std(Y(:));
stdz = std(Z(:));

D = [(Y'-muy)./stdy (Z'-muz)./stdz];
D = max(min(D,3.291),-3.291);

[~, grid] = kmeans(D, m, 'Start', 'uniform');
grid = [(stdy*grid(:,1)+muy) (stdz*grid(:,2)+muz)];
scatter(Y,Z); hold on; scatter(grid(:,1),grid(:,2),'filled'); hold off

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
YI     = Y(1:N)';
Y0     = [];
ZI     = Z(1:N)';
Z0     = []; 

while tol>10^-2 && count <= maxT

    % Simulate new observations
    Y      = zeros(N,nsim);
    Z      = zeros(N,nsim);
    epsy   = normrnd(0,1,[N nsim]);
    epsz   = normrnd(0,1,[N nsim]);

    Z(:,1) = ZI;
    Y(:,1) = YI;

    for ti=2:nsim
        Z(:,ti) = (1-rhoV)*eta    + rhoV.*Z(:,ti-1) + omega.*epsz(:,ti);
        Y(:,ti) = rhoM.*Y(:,ti-1) + sqrt(Z(:,ti)).*epsy(:,ti);
    end

    Y = real(Y'); Z=real(Z');
    [~,Gidx] = pdist2(grid,[Y(:) Z(:)],'euclidean','Smallest',1);
       Gidx  = reshape(Gidx,nsim,N)';
    Y = Y'; Z=Z';

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
    Y0     = [Y0 Y]; %#ok
    Z0     = [Z0 Z]; %#ok
    YI     = Y0(:,end);
    ZI     = Z0(:,end); 

    printtol = max(tol1,tol2);
    fprintf('Tolerance %4.6f Time T=%4.0f \n',printtol,count)
    tol      = printtol;

end

% Initial for EM 
[N,T]  = size(Y0);

sigY   = sqrt(mean((Y(:)-grid(Gidx,1)).^2,'all')).*ones(m,1);
sigZ   = sqrt(mean((Z(:)-grid(Gidx,2)).^2,'all')).*ones(m,1);

[~,Idx] = sort(grid(:,1)); 
gridy   = grid(Idx,1);
gridz   = grid(Idx,2);
Pi      = tpm1(Idx,Idx); 

pistat = ones(1,m)/(eye(m)-Pi+ones(m,m)); 
Y      = Y0';
Z      = Z0'; 
pis    = Pi;

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
    sigs      = sigY;
    phiy      = (exp(-((Y(1,:)-gridy).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    sigs      = sigZ;
    phiz      = (exp(-((Z(1,:)-gridz).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    phi{1}    = phiy.*phiz;

    alf{1}    = max(pistat'.*phi{1},10^-7);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        sigs       = sigY;
        phiy       = (exp(-((Y(ti,:)-gridy).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        sigs       = sigZ;
        phiz       = (exp(-((Z(ti,:)-gridz).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        phi{ti}    = phiy.*phiz;

        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},10^-7);
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

        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),10^-7);
        betscl{ti} = bet{ti}.*cscl{ti};

        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

    end
    ittoc = toc;
   
    %% Update Mean (i.e. Grid)
    munum   = zeros(m,1); 
    mudenom = zeros(m,1); 
    for ti=1:T
        munum      = munum + sum(gam{ti}.*Y(ti,:),2); 
        mudenom    = mudenom + sum(gam{ti},2);
    end 
    gridy  = munum./mudenom;

    munum   = zeros(m,1);
    mudenom = zeros(m,1);
    for ti=1:T
        munum      = munum + sum(gam{ti}.*Z(ti,:),2);
        mudenom    = mudenom + sum(gam{ti},2);
    end
    gridz = munum./mudenom;

    %% Update Transitions
    tpmnum  = zeros(m,m);
    for ti=1:(T-1)
        xis = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*pis;
        tpmnum  = tpmnum + xis;
    end 
    pis    = tpmnum./sum(tpmnum,2);
    pistat = ones(1,m)/(eye(m)-pis+ones(m,m)); % update initial

    %% Update sigY and sigZ  
    ssecumulY   = 0; 
    ssecumulZ   = 0; 
    for ti = 1:T
        ssetemp   = sum(gam{ti}.*((Y(ti,:)-gridy).^2),'all');
        ssecumulY = ssecumulY + ssetemp; 

        ssetemp   = sum(gam{ti}.*((Z(ti,:)-gridz).^2),'all');
        ssecumulZ = ssecumulZ + ssetemp; 
    end 
    sigY = sqrt(sum(ssecumulY,'all')./(N*T)).*ones(m,1);  
    sigZ = sqrt(sum(ssecumulZ,'all')./(N*T)).*ones(m,1); 

    % Check convergence 
    loglikelihood(niter,1) = ll;

    if niter > 1
        tol = abs(loglikelihood(niter,1)/loglikelihood(niter-1,1) -1);
    end

    fprintf('Iteration %4.0f Update Error %4.6f  Time %4.2f LogLik %4.6f \n',niter,tol,ittoc,ll)

    if tol<10^-3
        break
    end

end

%% output 
tpm   = pis; 
sigs  = [sigY(1) sigZ(1)]; 

end 

function [grid,tpm,ll,sigs,Y] = DiscretizeAR1(m,rho,sig)

%% Arguments
% m   - number of grid points
% rho - persistence
% sig - standard deviation of innovation

if m>49
    disp('Warning: Number of grid points exceeds 49. Fixing grid for numerical stability');
elseif mod(m,2)==0
    disp('Warning: Even number of grid points selected. Replacing with odd number m+1')
    m = m+1;
elseif abs(rho)>=1
    disp('ERROR: Absolute value of rho must be below 1')
    grid=[]; tpm =[]; ll=[]; sigs=[]; Y=[];
    return
end

%&------------------------------------------------------------------------%
%% Simulate Data & Initial Guess
%&------------------------------------------------------------------------%

% Grid
sigy   = sig/sqrt(1-rho^2);
space  = min(4,1.2*log(m-1));
grid   = linspace(-space*sigy,space*sigy,m)';

% Transition
Pi = zeros(m,m);
for ri=1:m
    Pi(ri,:) = normpdf(grid,rho*grid(ri),sig);
end
Pi = Pi./sum(Pi,2);

%_________________________________________________________________________%
% Recursive binning

N    = 100;
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
YI     = normrnd(0,sigy,[N 1]);
Y0     = [];


%_________________________________________________________________________%
% Update transitions recursively

while tol>10^-2 && count <= maxT

    % Simulate new observations
    eps    = normrnd(0,sig,[N nsim]);
    Y      = zeros(N,nsim);
    Y(:,1) = rho.*YI + eps(:,1);
    for ti = 2:nsim
        Y(:,ti) = rho.*Y(:,ti-1) + eps(:,ti);
        Y(:,ti) = min(max(Y(:,ti),-4*sigy),4*sigy);
    end

    Y = Y';
    [~,Yidx] = pdist2(grid,Y(:),'euclidean','Smallest',1);
    Yidx     = reshape(Yidx,nsim,N)';
    Y = Y';

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
    E    = flip(eye(m));
    C    = E*tpm1*E;
    tpm1 = (tpm1+C)./2;
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
pistat    = ones(1,m)/(eye(m)-Pi+ones(m,m));
sigs      = ones(m,1).*sqrt(mean((Y-grid(Yidx)).^2,'all'));
Y         = Y0';
pis       = Pi;

%&------------------------------------------------------------------------%
%%  EM Algorithm
%&------------------------------------------------------------------------%

% Forward and backward probabilities
alf     = cell(T,1);
alfscl  = cell(T,1);
bet     = cell(T,1);
betscl  = cell(T,1);
phi     = cell(T,1);
gam     = cell(T,1);
cscl    = cell(T,1);

maxiter       = 1000;
loglikelihood = zeros(maxiter,1);

inititer = 1; % initialize the sigma before updating the mean and variance

for niter = 1:maxiter

    % forward
    phi{1}    = (exp(-((Y(1,:)-grid).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
    alf{1}    = max(pistat'.*phi{1},1e-8);
    cscl{1}   = 1./sum(alf{1},1);
    alfscl{1} = alf{1}.*cscl{1};
    llk       = -log(cscl{1});
    ll        = sum(llk,"all");

    for ti = 2:T

        phi{ti}    = (exp(-((Y(ti,:)-grid).^2)./(2*sigs.^2)))./(sqrt(2*pi).*sigs);
        alf{ti}    = max(((pis')*alfscl{ti-1}).*phi{ti},1e-8);
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

        bet{ti}    = max(pis*(betscl{ti+1}.*phi{ti+1}),1e-8);
        betscl{ti} = bet{ti}.*cscl{ti};
        gtemp      = betscl{ti}.*alfscl{ti};
        gam{ti}    = gtemp./sum(gtemp,1);

    end

    % Update Mean (i.e. Grid), imposing symmetry
    if m<50 && niter > inititer
        munum   = zeros(m,1);
        mudenom = zeros(m,1);
        for ti=1:T
            munum      = munum + sum(gam{ti}.*Y(ti,:),2);
            mudenom    = mudenom + sum(gam{ti},2);
        end
        grid  = sort((munum - flip(munum))./(mudenom+flip(mudenom)));
    end

    % Update Transitions, imposing symmetry
    if niter > inititer
        tpmnum  = zeros(m,m);
        for ti=1:(T-1)
            xis = alfscl{ti}*(phi{ti+1}.*betscl{ti+1})'.*pis;
            tpmnum  = tpmnum + xis;
        end
        Pi     = tpmnum;
        E      = flip(eye(m));
        C      =  E*Pi*E;
        Pi     = (Pi+C)./2;
        pis    = Pi./sum(Pi,2);
    end
    pistat = max(ones(1,m)/(eye(m)-pis+ones(m,m)),1e-7); % update initial

    % Update sigma
    ssecumul = 0;
    for ti = 1:T
        ssetemp  = sum(gam{ti}.*((Y(ti,:)-grid).^2),'all');
        ssecumul = ssecumul + ssetemp;
    end
    sigs = sqrt(sum(ssecumul,'all')./(N*T)).*ones(m,1);

    % Check convergence
    loglikelihood(niter,1) = ll;
    if niter > inititer+1
        tol = abs(ll/loglikelihood(niter-1,1) -1);
        fprintf('Iteration %4.0f Update Error %4.6f LogLik %4.6f \n',niter-inititer,tol,ll/(N*T))
        if tol< 0.005
            break
        end
    end

end

%% output
tpm  = pis;
sigs = sigs(1);
ll   = loglikelihood(niter)/(N*T);

end

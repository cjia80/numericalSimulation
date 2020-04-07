
% support recovery probability for   LQ-LS(PDASC) 
% on different sparsity level
clc
clear all 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
n = 1000;        % Signal dimension
m = 500;         % Number of measurements
rflip = 0.01;    % probability of sign flips
sigma = 0.05;
rho = 0.1;
KK = 1:2:20;
len = length(KK);
numrun = 100;
total_err = zeros(len,numrun);
total_time = zeros(len,numrun);
total_supp = zeros(len,numrun);
for k = 1:len
    k
    K  = KK(k);
    for ii = 1 : numrun
        ii   
        % --------- Generate data ---------------------
        seednum = k*ii+1000;
        rand('seed',seednum);   % fix seed 
        randn('seed', seednum); % fix seed 
        xtrue = zeros(n,1);
        rp = randperm(n);
        xtrue(rp(1:K)) = sign(randn(K,1));
        xtrue = xtrue/norm(xtrue);
        supp = find(xtrue);
        SIGMA = zeros(n,n);
        for kk = 1:n
            for jj = 1:n
               SIGMA(kk,jj) = rho^(abs(kk-jj));
            end
        end
        Mu = zeros(1,n);
        Psi =  mvnrnd(Mu,SIGMA,m);
        %Psi = randn(m,n);
        ye = Psi*xtrue;
        noise = sigma*randn(m,1);
        y = sign( ye + noise);
        nflip = floor(rflip*m);
        indxflip = randperm(m);
        indxflip = indxflip(1:nflip);
        y(indxflip) = -y(indxflip);
    %--------------------------------------------------------
        fprintf('------ LQ-LS(wPDASC) ------\n')
        tic
           [x,lam,ithist] = wpdasc(Psi,Psi',y);
           esupp = find(x);
           Psiesupp = Psi(:,esupp); 
           % debias
           x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
           total_time(k,ii) =  toc;
           total_err(k,ii) =  norm(xtrue - x);
           rs = setdiff(esupp,supp);
           rse = setdiff(supp,esupp);
           if isempty(rs) && isempty(rse)
                total_supp(k,ii) = 1;
                
           end
    end
end
% show results 
avsupp = mean(total_supp,2);
averror = mean(total_err,2);

figure(1)
% subplot(1,2,1)
plot(KK,avsupp,'ro--','LineWidth',1.5)
h = xlabel('$S$');
set(h,'Interpreter','latex','fontsize',15)
set(gca,'xtick',0:2:20)
h = ylabel(' probability ');
set(gca,'ytick',0:0.2:1)
set(h,'Interpreter','latex','fontsize',15)
axis([0 20 0 1.1])

% subplot(1,2,2)
plot(KK,averror,'ro--','LineWidth',1.5)
h = xlabel('$S$');
set(h,'Interpreter','latex','fontsize',15)
set(gca,'xtick',0:2:20)
h = ylabel(' $\ell_2$-error ');
set(gca,'ytick',0:0.1:0.6)
set(h,'Interpreter','latex','fontsize',15)
axis([0 20 0 0.7])
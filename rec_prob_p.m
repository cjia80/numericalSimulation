clc
clear all 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
n = 1000;        % Signal dimension
m = 500;         % Number of measurements
K = 5;           % signal sparsity
rho = 0.1;       % coorelation
sigma = 0.05;
 Rflip = 0:0.03:0.15;
% Rflip = 0.85:0.03:1;
len = length(Rflip);
numrun = 100;
total_err = zeros(len,numrun);
total_time = zeros(len,numrun);
total_supp = zeros(len,numrun);
for k = 1:len
    k
    rflip = Rflip(k);    
    for ii = 1 : 100
        ii   
        % --------- Generate data ---------------------
        seednum = k*ii;
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
        ye = Psi*xtrue;
        noise = sigma*randn(m,1);
        y = sign( ye + noise);
        nflip = floor(rflip*m);
        indxflip = randperm(m);
        indxflip = indxflip(1:nflip);
        y(indxflip) = -y(indxflip);
    %--------------------------------------------------------
        fprintf('------ LQ-LS(WPDASC) ------\n')
        tic
        [x,lam,ithist] = wpdasc(Psi,Psi',y);
        toc
           % debias
           esupp = find(x);
           Psiesupp = Psi(:,esupp); 
           x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
           total_time(k,ii) =  toc;
           total_err(k,ii) =  norm(xtrue - x);
           esupp = find(x);
            if length(esupp) == length(supp)
                rs = setdiff(esupp,supp);
                if isempty(rs)
                     total_supp(k,ii) = 1;
                end
            end
    end
end
% show results 
avsupp = mean(total_supp,2);
averror = mean(total_err,2);

figure(1)
% subplot(1,2,1)
plot(Rflip,avsupp,'ro--','LineWidth',1.5)
h = xlabel('$p$');
set(h,'Interpreter','latex','fontsize',15)
set(gca,'xtick',0:0.03:0.15)
h = ylabel(' Probability ');
set(gca,'ytick',0:0.2:1)
set(h,'Interpreter','latex','fontsize',15)
axis([0 0.15 0 1.2])

% subplot(1,2,2)
% plot(Rflip,averror,'ro--','LineWidth',1.5)
% h = xlabel('$p$');
% set(h,'Interpreter','latex','fontsize',15)
% set(gca,'xtick',0:0.03:0.15)
% h = ylabel(' $\ell_2$-error ');
% set(gca,'ytick',0:0.05:0.25)
% set(h,'Interpreter','latex','fontsize',15)
% axis([0 0.15 0 0.25])
% 
% % 
% % subplot(1,2,1)
% % plot(Rflip,avsupp,'ro--','LineWidth',1.5)
% % h = xlabel('$p$');
% % set(h,'Interpreter','latex','fontsize',15)
% % set(gca,'xtick',0.85:0.03:1)
% % h = ylabel(' Probability ');
% % set(gca,'ytick',0:0.2:1)
% % set(h,'Interpreter','latex','fontsize',15)
% % axis([0.85 1 0 1.2])
% 
% subplot(1,2,2)
% plot(Rflip,averror,'ro--','LineWidth',1.5)
% h = xlabel('$p$');
% set(h,'Interpreter','latex','fontsize',15)
% set(gca,'xtick',0.85:0.03:1)
% h = ylabel(' $\ell_2$-error ');
% set(gca,'ytick',1.95:0.01:2)
% set(h,'Interpreter','latex','fontsize',15)
% axis([0.85 1 1.95 2])


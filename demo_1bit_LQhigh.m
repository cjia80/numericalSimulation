clc
clear all 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data 
seednum = 100;
rand('seed',seednum);   % fix seed 
randn('seed', seednum); % fix seed 
m = 500;                % number of samples
n = 1000;               % number of dimention
K =5;                  % sparsity level
sigma = 0.1;           % noise level
nflip = 10;             % number of sign fips
xtrue = zeros(n,1);
q = randperm(n);
supp = q(1:K); 
xtrue(supp) = sign(randn(K,1));  
xtrue = xtrue/norm(xtrue);
supp = find(xtrue);
rho = 0.1;
SIGMA = zeros(n,n);
for k = 1:n
    for j = 1:n
       SIGMA(k,j) = rho^(abs(k-j));
    end
end
Mu = zeros(1,n);
Psi =  mvnrnd(Mu,SIGMA,m);
noise = randn(m,1);
ye = Psi*xtrue + sigma*noise;
xtrue = xtrue/norm(xtrue);
y = sign(ye);
id = randperm(m);
idflip = id(1:nflip);
yidflip  = y(idflip);
y(idflip) = - yidflip;
% Phi = Psi;
%  Phit = Psi';
%  Psi =  @(z) Phi*z;
% Psit =  @(z) Phit*z;
implicit = isa(Psi,'function_handle');

if implicit == 0
    tic
   [xb,lamb,ithistb] = wpdasc(Psi,Psi',y);
  
   toc 
    
     % debias
      esuppb=find(xb);
       Psiesuppb = Psi(:,esuppb); 
       xb(esuppb) = (Psiesuppb'*Psiesuppb)\(Psiesuppb'*y);
       xb = xb/norm(xb);
       [a,b]=sort(abs(xtrue),'descend');
       aa=a(1:K);
       [a1,b1]=sort(abs(xb),'descend');
       aa1=a1(1:K);
    figure(1) 
    plot(1:n,xtrue,'ko',1:n,xb,'r*'), 
     h = xlabel(' the $i$-coordinate of a vector');
set(h,'Interpreter','latex','fontsize',12)
h = ylabel('the $i$-th value of a vector');
set(h,'Interpreter','latex','fontsize',12)
h=legend({'$x^*$','$x_{\hat{\lambda}}$'},'Interpreter','latex','fontsize',10);
set(h,'Interpreter','latex','fontsize',10)
axis([0 1000 -0.5 0.5])
    
     figure(2),
    plot(ithistb.as,'r+','LineWidth',2), 
    h=xlabel('the index $k$ of $\lambda_k$ ');
    set(h,'Interpreter','latex','fontsize',12)
    h=ylabel({'$\|x_{\lambda_k}\|_{0}$'},'Interpreter','latex','fontsize',12);
set(h,'Interpreter','latex','fontsize',12)

else
      tic
       [xb,lamb,ithistb] = wpdasc(Psi,Psit,y);
       xb = xb/norm(xb);
       % debias
       esuppb = find(xb);
       rhs = Psit(y);
       option.subset = esuppb;
       option.n = n;
       option.init = xb(esuppb);
       option.itcg = length(esuppb);
       xb(esuppb) = Subcg(Psi,Psit,rhs(esuppb), option);
       xb = xb/norm(xb);
      toc
      
     figure(1) 
    plot(1:n,xtrue,'ko',1:n,xb,'r*'), 
     h = xlabel(' the $i$-th of a vector');
    set(h,'Interpreter','latex','fontsize',14)
    h = ylabel('the $i$-th element of a vector');
    set(h,'Interpreter','latex','fontsize',14)
    h=legend({'$x^*$','$x_{\hat{\lambda}}$'},'Interpreter','latex','fontsize',14);
   set(h,'Interpreter','latex','fontsize',14)
   axis([0 1000 -0.5 0.5])

    
     figure(2),
    plot(ithistb.as,'r+','LineWidth',2), 
    h=xlabel('the nonzero active set path');
    set(h,'Interpreter','latex','fontsize',14)
    h=ylabel({'$x_{\hat{\lambda}}$'},'Interpreter','latex','fontsize',14);
    set(h,'Interpreter','latex','fontsize',14)
end
fdb = setdiff(esuppb,supp); % false dis
fmb = setdiff(supp,esuppb); % false mis
err_xb = norm(xtrue - xb)

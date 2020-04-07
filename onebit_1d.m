clc
clear all 
close all
addpath(genpath(fileparts(mfilename('fullpath'))));
warning off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate data 
seednum = 0;
signalsize = 8000;
nsam = 2500; 
sigma = 5e-1;
[Psi,Psit,y,xtrue,supp,out] = getrealdata(signalsize,nsam,sigma,seednum);
m = length(y);      % number of samples
n = length(xtrue);  % number of dimention
K = length(supp);   % sparsity level
nflip = 100;        % number of sign fips
id = randperm(m);
idflip = id(1:nflip);
yidflip  = y(idflip);
y(idflip) = - yidflip;
xtrue = xtrue/norm(xtrue);
scale = max(max(out.I)) - min(min(out.I));
mm=8;
maxiter=100;
totaltime=zeros(maxiter,mm);
totalerr=zeros(maxiter,mm);
totalpsnr=zeros(maxiter,mm);


fprintf('------ BIHT ----------\n')
tic
x = BIHT_1(y, Psi, K,Psit);
bihttime = toc
err_biht = norm(x - xtrue)
xbiht = out.W(x) + out.bar;
figure(2)
plot(out.mesh,xbiht,'k-','LineWidth',1)
box off
axis off  
psnrbiht = psnr(out.I,xbiht,scale);
title(['BIHT, PSNR = ',num2str(psnrbiht,'%2.0f')])


fprintf('------ AOP  ------\n')
L  = nflip;  % number of wrong labels.
alpha   = 1;
tic
    x = BIHT_AOP_flip(y, Psi, Psit(y), K, L, 1, 100, alpha,Psit);
    x = x/norm(x);
aoptime = toc
err_aop = norm(x - xtrue)
xaop = out.W(x) + out.bar;
figure(3)
plot(out.mesh,xaop,'k-','LineWidth',1)
box off
axis off  
psnraop = psnr(out.I,xaop,scale);
title(['AOP, PSNR = ',num2str(psnraop,'%2.0f')])


fprintf('------ LinProj  ------\n')
tic 
    tx = Psit(y)/m;
    x = LinProj(tx,sqrt(K),1,.5,100,50,1e-6);
linprojtime = toc
err_linproj = norm(x - xtrue)
xlinproj = out.W(x) + out.bar;
figure(4)
plot(out.mesh,xlinproj,'k-','LineWidth',1)
box off
axis off  
psnrlinproj = psnr(out.I,xlinproj,scale);
title(['LP, PSNR = ',num2str(psnrlinproj,'%2.0f')])



fprintf('------ PBAOP  ------\n')
L = nflip;      % using  number sign flips.
alpha   = 1;
tau = 0.05;
tic
    x = PIHT_AOP_flip(y, Psi, Psit(y), K, L, 1, 100, alpha, tau,Psit);
    x = x/norm(x);
paoptime = toc
err_paop = norm(x - xtrue)
xpaop = out.W(x) + out.bar;
figure(5)
plot(out.mesh,xpaop,'k-','LineWidth',1)
box off
axis off  
psnrpaop = psnr(out.I,xpaop,scale);
title(['PBAOP, PSNR = ',num2str(psnrpaop,'%2.0f')])


fprintf('------ PDASC  ----------\n')
tic
   [x,lam,ithist] = pdasc(Psi,Psit,y);
   x = x/norm(x);
hopdastime = toc
err_hopdas = norm(x - xtrue)
xhopdas = out.W(x) + out.bar;
figure(6)
plot(out.mesh,xhopdas,'k-','LineWidth',1)
box off
axis off  
psnrhopdas = psnr(out.I,xhopdas,scale);
title(['PDASC, PSNR = ',num2str(psnrhopdas,'%2.0f')])

fprintf('------ WPDASC ----------\n')
  tic
   [x,lam,ithist] = wpdasc(Psi,Psit,y);
   x = x/norm(x);
wpdastimeb = toc
err_wpdasb = norm(x - xtrue)
wpdas = out.W(x) + out.bar;
figure(7)
plot(out.mesh,wpdas,'k-','LineWidth',1)
box off
axis off  
psnrhopdas = psnr(out.I,wpdas,scale);
title(['WPDASC, PSNR = ',num2str(psnrhopdas,'%2.0f')])

 
  fprintf('------ Passive ----------\n')
  tic
    x=passive_1bit(y,Psi,xtrue,K,Psit);
   x = x/norm(x);
passtimeb = toc
err_passiveb = norm(x - xtrue)
xhopass = out.W(x) + out.bar;
figure(8)
plot(out.mesh,xhopass,'k-','LineWidth',1)
box off
axis off  
psnrhopass = psnr(out.I,xhopass,scale);
title(['Passive, PSNR = ',num2str(psnrhopass,'%2.0f')])


fprintf('------ MCP ----------\n')
  tic
 x = mcp_1bit(y,Psi,xtrue,K,Psit);
   x = x/norm(x);
mcptimeb = toc
err_mcpb = norm(x - xtrue)
xhomcp = out.W(x) + out.bar;
figure(9)
plot(out.mesh,xhomcp,'k-','LineWidth',1)
box off
axis off  
psnrhomcp = psnr(out.I,xhomcp,scale);
title(['MCP, PSNR = ',num2str(psnrhomcp,'%2.0f')])

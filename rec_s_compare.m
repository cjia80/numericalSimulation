clc
clear all 
close all
warning off
addpath(genpath(fileparts(mfilename('fullpath'))));
% All available solvers:
FIELDS={'BIHT','AOP','LinProj','PBAOP','PDASC','WPDASC','Passive','MCP'};
n = 1000;        % Signal dimension
m = 500;         % Number of measurements
rflip = 0.01;    % probability of sign flips
sigma = 0.05;
rho = 0.1;
KK =1:2:16;
% sigma=0.5;   
% rflip=0.05;
% rho=0.5;
% KK =1:2:16;
len = length(KK);
maxnumtest = 100;
nmethod=8;
total_time = zeros(maxnumtest,len,nmethod);     
total_err = zeros(maxnumtest,len,nmethod);
total_supp= zeros(maxnumtest,len,nmethod); 

for k = 1:len
    k
    K  = KK(k)
   for ii = 1 : maxnumtest
        ii   
%         --------- Generate data ---------------------
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
%         implicit = isa(Psi,'function_handle');

        %Psi = randn(m,n);
        ye = Psi*xtrue;
        noise = sigma*randn(m,1);
        y = sign( ye + noise);
        nflip = floor(rflip*m);
        indxflip = randperm(m);
        indxflip = indxflip(1:nflip);
        y(indxflip) = -y(indxflip);
  %---------------------BIHT------------------------------------   
    ff='BIHT';
   if ismember(ff,FIELDS) 
    fprintf('------BIHT ------\n')
    tic
    x = BIHT_1(y, Psi, K);
    bihttime=toc;
    mm=1;    
    total_time(ii,k,mm) = bihttime;
    total_err(ii,k,mm) =  norm(xtrue - x);
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,k,mm) =  1;
        end
    end
   end
   %---------------------AOP------------------------------------   
   ff='AOP';
   if ismember(ff,FIELDS) 
    fprintf('------AOP ------\n')
    L  = nflip;  % using the real number of wrong labels.
    alpha   = 1;
    tic
    x = BIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha);
    x = x/norm(x);
    aoptime=toc;
    mm=2;
    total_time(ii,k,mm) = aoptime;
    total_err(ii,k,mm) =  norm(xtrue - x);
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,k,mm) = 1;
        end
    end
   end  
    %-------------------------------------------------------
    ff='LinProj';
    if ismember(ff,FIELDS) 
    fprintf('------ LinProj ------\n')
       tic,
     tx = Psi'*y/m;
     x = LinProj(tx,sqrt(K),1,.5,100,50,1e-6);
     lptime=toc;
     mm=3;
     total_time(ii,k,mm) =  lptime;
     total_err(ii,k,mm) =  norm(xtrue - x);    
     esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,k,mm) = 1;
        end
    end
    end
    %---------------------------------------------------
  ff='PBAOP'; 
   if ismember(ff,FIELDS) 
    fprintf('------ PBAOP ------\n')
     L = nflip;      % using  number of sign filips.
    alpha   = 1;
    tau = 0.05;
    tic
    x = PIHT_AOP_flip(y, Psi, Psi'*y, K, L, 1, 100, alpha, tau);
    x = x/norm(x);
    pbtime=toc;
    mm=4;
   total_time(ii,k,mm) = pbtime;
    total_err(ii,k,mm) =  norm(xtrue - x);
    esupp = find(x);
    if length(esupp) == length(supp)
        rs = setdiff(esupp,supp);
        if isempty(rs)
             total_supp(ii,k,mm) = 1;
        end
    end
   end
     %--------------------------------------------------------
  ff='PDASC';
   if ismember(ff,FIELDS) 
         fprintf('------L1-LS(PDASC) ------\n')
        tic
           [x,lam,ithist] = pdasc(Psi,Psi',y);
           esupp = find(x);
           Psiesupp = Psi(:,esupp); 
           % debias
           x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
           x = x/norm(x);
           pdastime=toc;
           mm=5;
           total_time(ii,k,mm) =  pdastime;          
           total_err(ii,k,mm) =  norm(xtrue - x);
           rs = setdiff(esupp,supp);
           rse = setdiff(supp,esupp);
           if isempty(rs) && isempty(rse)
                total_supp(ii,k,mm) = 1;
           end
   end
     
    
 % ------------------------------------------------
    ff='Passive';
      if ismember(ff,FIELDS) 
       fprintf('------ Passive------\n')
        mu = sqrt(log(n)/m);
       tic
       x=passive_1bit(y,Psi,xtrue,K); 
        passtime=toc;
        mm=6;
         total_time(ii,k,mm) =  passtime;          
         total_err(ii,k,mm) =  norm(xtrue - x);
         esupp = find(x);
         if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
             if isempty(rs)
             total_supp(ii,k,mm) = 1;
             end
         end
      end
      
    % ------------------------------------------------
    ff='MCP';
      if ismember(ff,FIELDS) 
       fprintf('------ MCP------\n')
        mu = sqrt(log(n)/m);
       tic
       x=mcp_1bit(y,Psi,xtrue,K); 
        passtime=toc;
        mm=7;
         total_time(ii,k,mm) =  passtime;          
         total_err(ii,k,mm)=  norm(xtrue - x);
         esupp = find(x);
         if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
             if isempty(rs)
             total_supp(ii,k,mm) = 1;
             end
         end
      end   
   %--------------------------------------------------------
   ff='WPDASC';
   if ismember(ff,FIELDS) 
   fprintf('------ WPDASC------\n')
         tic
        [x,lam,ithist] = wpdasc(Psi,Psi',y);
        esupp = find(x);
         Psiesupp = Psi(:,esupp); 
         % debias
         x(esupp) = (Psiesupp'*Psiesupp)\(Psiesupp'*y);
         x = x/norm(x);
         pdasbtime=toc;
         mm=8;
         total_time(ii,k,mm) =  pdasbtime;          
         total_err(ii,k,mm) =  norm(xtrue - x);
         esupp = find(x);
       if length(esupp) == length(supp)
             rs = setdiff(esupp,supp);
          if isempty(rs)
             total_supp(ii,k,mm) = 1;
          end
        end
   end    
   end
end
averagetime = mean(total_time);
averageerror = mean(total_err);
proboracle = mean(total_supp);
K = KK;

prob = [proboracle(:,:,1);proboracle(:,:,2);proboracle(:,:,3);proboracle(:,:,4);proboracle(:,:,5);proboracle(:,:,6);proboracle(:,:,7);proboracle(:,:,8)]';%proboracle(:,:,9)
figure(1)
plot(K,prob(:,1),'b*:',K,prob(:,2),'ro--',K,prob(:,3),'x-.',K,prob(:,4),'gs:',K,prob(:,5),'kd-',K,prob(:,6),'bp-',K,prob(:,7),'ch-.',K,prob(:,8),'m+-','LineWidth',2)
h = xlabel('$S$');
set(h,'Interpreter','latex','fontsize',12)
set(gca,'xtick',0:2:16)
h = ylabel('Probability');
set(gca,'ytick',0:0.2:1.4)
set(h,'Interpreter','latex','fontsize',12)
h = legend({'BIHT','AOP','LinProj','PBAOP','PDASC','Passive','MCP','WPDASC'},'FontSize',8);
set(h,'Interpreter','latex','fontsize',8)
axis([0 16 0 1.6])



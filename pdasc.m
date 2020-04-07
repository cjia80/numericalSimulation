function [x,lam,ithist] = pdasc(X,Xt,y,opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             min      1/(2m)||X*x-y||^2  + lambda ||x||_1              %
%    by PDAS  algorithm with  continuation Lam ={lam_{1},...,lam_{N}    %
%              (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)           %
%    Created on Oct 17, 2016                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                      %
% INPUTS:                                                               %
%     X  ---  sensing matrix (R^{m*n})                                  %
%     Xt ---  transpose of X                                            %
%     y  ---  data vector                                               %
%   opts ---  structure containing                                      %
%         N   --- length of path (default: 200)                         %
%       Lmin  --- minimun in Lam (default: 1e-4)                        %
%        mu   --- stop if \|x_{lam_k}\|_0 > mu                          %
%                 (default:min(0.5*m/log(n),sqrt(n)))                   %
%       init  --- initial value for  (default: 0)                       %   
% OUTPUTS:                                                              %
%     x   ---- recovered signal                                         %
%    lam  ---- regularization  parameter                                %
%  ithist ---- structure on iteration history, containing               %
%          .x  --- solution path                                        %
%          .as --- size of active set  on the path                      %
%          .it --- # of iteration on the path                           %
% ======================================================================%
 implicit = isa(X,'function_handle');
if implicit == 0
   [m,n] = size(X);
   Xty = Xt*y;
else
   m = length(y);
   Xty = Xt(y);
   n = length(Xty);
end
linf  = norm(Xty,inf);
if ~exist('opts','var')
    opts.N     = 200;
    opts.mu    = min(0.5*m/log(n),sqrt(n));
    opts.Lmax = 1;
    opts.Lmin  = 1e-4;
    opts.init  = zeros(n,1);
    opts.n = n;
    opts.m = m;
    opts.Xty  = Xty;
end
% construct the homotopy path
Lam   = exp(linspace(log(opts.Lmax),log(opts.Lmin),opts.N))';
Lam   = Lam(2:end);
Lam = Lam*linf/m;
ithist.Lam = Lam;
%% main loop for pathfolling and choosing lambda and output solution
ithist.x   = [];
ithist.as  = [];
for k = 1:length(Lam)
    opts.lam = Lam(k);
    [x,s] = pdas(X,Xt,opts);
    opts.init = x;
    ithist.x(:,k) = x;
    ithist.as = [ithist.as;s];    % size of active set     
    if s > opts.mu  
          % display('# NON-ZERO IS TOO MUCH, STOP ...')
          break   
    end
end
% select the solution on the path by voting 
ii = find(ithist.as == mode(ithist.as));
ii = ii(end);
x = ithist.x(:,ii);
lam = Lam(ii);
end
%% sub functions
function [x,s] = pdas(X,Xt,opts)
%-------------------------------------------------------------------------%
%         Solving                                                         %
%                    1/(2m)||X*x-y||^2  + lambda ||x||_1                  %
%         by  one step primal-dual active set algorithm                   %
%-------------------------------------------------------------------------%
% INPUTS:                                                                 %
%     X  ---  sensing matrix (R^{m*n})                                    %
%     Xt ---  transpose of X                                              %
%     y  ---  data vector                                                 %
%   opts ---  structure containing                                        %
%          lam  ---  regualrization paramter                              %
%         init  ---  initial guess                                        %
% OUTPUTS:                                                                %
%      x  ---  solution                                                   %                                                    
%      s  ---  size of active set                                         %
%              (c) by Yuling Jiao (yulingjiaomath@whu.edu.cn)             %
%    Created on Oct 17, 2016                                              %
%-------------------------------------------------------------------------%
lam = opts.lam;
Xty = opts.Xty;
x0 = opts.init;
n  = opts.n;
m  = opts.m;
implicit = isa(X,'function_handle');
if implicit == 0
      % initializing ...
      pd  = x0 + (Xty - Xt*(X*x0))/m;
      A  = find(abs(pd)>lam);
      s = length(A);
      x = zeros(n,1);
      dA = lam*sign(pd(A));
      XtyA = Xty(A);
      rhs = XtyA - m*dA;
      XA = X(:,A);
      G = XA'*XA;
      x(A) = G\rhs;
else
   % initializing ...
    pd  = x0 + (Xty - Xt(X(x0)))/m;
   A  = find(abs(pd)>lam);
   s = length(A);
   x = zeros(n,1);
   dA = lam*sign(pd(A));
   XtyA = Xty(A);
   rhs = XtyA - m*dA;
   option.subset = A;
   option.n = n;
   option.init = x0(A);
   option.itcg = min(s,2); 
   x(A) = Subcg(X,Xt, rhs,option);
end
end

function  x = Subcg(L,Lt,b, option)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve a symmetric positive definite system (Lt*L)x = b via CG           %
% Input:                                                                 %
%            L - Either an    matrix, or a function handle               %
%            Lt- Transpose of L                                          %
%            b - vector                                                  %
%      maxiter - Maximum number of iterations (defaut length(b))         %
%      initial - Initial value (defaut 0)                                %
%       option - a stucture for subset cg                                %
%             .subset    -  index of conlums of A used                   %
%             .n         -  number of whole conlums                      %
%             .init      -  initial value                                %
%             .itcg      -  number of iterations                         %
% Output:                                                                %
%            x -  soultion                                               %
% Copyright (c)  Yuling Jiao(yulingjiaomath@whu.edu.cn)                  %
% Created on 17 October, 2013                                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = option.n;
maxiter = option.itcg;
subset = option.subset;
Ssubt = @(z) upsam(z,subset,n);
Ssub = @(z) z(subset,:);
Aop = @(z) Ssub(Lt(L(Ssubt(z))));
x = option.init; 
r = b - Aop(x); 
d = r;
delta = r'*r;
iter = 0;
while iter < maxiter
  q = Aop(d); 
  Alpha = delta/(d'*q);
  x = x + Alpha*d;
  r = r - Alpha*q;
  deltaold = delta;
  delta = r'*r;
  beta = delta/deltaold;
  d = r + beta*d;
  iter = iter + 1;
end
end
%% subfunctions
function upz = upsam(z,id,nn)
  upz = zeros(nn,1);
  upz(id) = z;
end
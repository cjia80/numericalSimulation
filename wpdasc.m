function [x,lam,ithist] = wpdasc(X,Xt,y,opts)
%-----------------------------------------------------------------------%
% Minimizing  1/2||X*x-y||^2 + lambda ||x||_q^q                           %
% by PDASC algorithm, with automatically choosing lambda                %
% INPUTS:                                                               %
%     X  ---  sampling matrix (normalized)                              %
%     Xt ---  transpose of sampling matrix                              %
%     y  ---  data vector                                               %
%   opts ---  structure containing                                      %
%         N   --- length of path (defaut: 100)                          %
%       Lmin  --- minimun in Lam (defaut: 1e-10)                        %
%        del  --- noise level                                           %
%        mu   --- stop if ||x_{lambda_{k}}||> opts.mu*n (defaut: 0.5)   %
%        x0   --- initial value for PDASC (defaut: 0)                   %
%     MaxIter --- maximum number of iterations in PDAS  (defaut: 5)     %
% OUTPUTS:                                                              %
%     x   ---- reconstructed signal                                     %
%    lam  ---- regu. parameter (by discrepancy principle or BIC)        %
%  ithist ---- structure on iteration history, containing               %
%          .x  --- solution path                                        %
%          .as --- size of active set                                   %
%          .it --- # of iteration on the path                           %
%          .res--- residual on the path                                 %
%------------------------------------------------------------------------%
% Based on huang2013 UPDASC                            %
% Created on Jan, 2020                                               %
% ======================================================================%
implicit = isa(X,'function_handle');
if implicit == 0
    [n,p] = size(X);
    Xty = Xt*y;
else
    n = length(y);
    Xty = Xt(y);
    p= length(Xty);
end

if ~exist('opts','var')
    opts.N     = 200;
    opts.Lmax  = 1;
    opts.Lmin  = 1e-4;
    opts.MaxIt = 1;
    opts.mu    = min(0.5*n/log(n),sqrt(n));
    opts.n = n;
    opts.p = p;
    opts.Xty  = Xty;
     opts.weight = 0.6;
      opts.tau = 0.9;   
      opts.del = 1.5*sqrt(n)*0.5;
      opts.sel = 'vote';   % you can set it to 'bic'
end
tau=  opts.tau;
linf  = norm(Xty,inf);  
cnst = (linf/(2-tau)*(2*(1-tau))^((tau-1)/(2-tau)))^(2-tau);
Lam  = exp(linspace(log(opts.Lmax),log(opts.Lmin),opts.N))';
Lam  = Lam(2:end);
Lam  = Lam*cnst/n;
ithist.Lam = Lam;
%% main loop for pathfolling and choosing lambda and output solution
ithist.x   = [];
ithist.as  = [];
ithist.it  = [];
ithist.res = [];
ithist.bic = [];
ny = norm(y)^2;
opts.Xty = Xty;
opts.G  = [];
opts.A = [];
d0 =Xty/n;   % initial with 0
x0 = zeros(p,1);
for k = 1:length(Lam)
     opts.lam = Lam(k);
   [x,d,G,A,s,L,It] = pdas(X,Xt,y,opts,x0,d0);
    % warmstarts 
    x0 = x;
    d0 = d;      
    opts.G = G;
    opts.A = A;
    ithist.x = [ithist.x, x];
 ithist.it = [ithist.it;It];  % iter #
    ithist.as = [ithist.as;s];    % size of active set
    res       = ny - L;                
    ithist.res= [ithist.res; res];
    ithist.bic= [ithist.bic; log(res/n) + log(n)*s/n];
    if res <= opts.del
        % discrepancy principle
        % disp('Discrepancy principle is satisfied.')
        lam = Lam(k);
        x   = ithist.x(:,k);
        id = k;
        break;
    end
    if s > opts.mu
        % display('Exceed maximum degree of freedoms')
        lam = Lam(k);
        break
    end
end
if res > opts.del
   switch opts.sel
       case 'vote'  % Voting
            % display('Selection by voting')
            id = find(ithist.as==0); % disimiss the case that the active set=0
            ithist.as(id) = [];
            ithist.x(:,id)=[];
            Lam(id) = [];
            ii = find(ithist.as == mode(ithist.as));
            id = ii(end);
            x  = ithist.x(:,id);
            lam= Lam(id);
       otherwise      % Bayesian information criterion
            % display('Selection by BIC')
            id = find(ithist.as==0);
            ithist.as(id) = [];
            ithist.x(:,id)=[];
            Lam(id) = [];
            ithist.bic (id) = [];
            ii = find(ithist.bic == min(ithist.bic));
            id = ii(end);
            x  = ithist.x(:,id);
            lam= Lam(id);
    end
end
end
% subfunction
function [x,d,G,A,s,L,It] = pdas(X,Xt,y,opts,x0,d0)
%-------------------------------------------------------------------------%
% Solving nonconvex problem                                               %
%      1/2||X*x-y||^2 + _\lam\|x\|^q_q                 %
% by a primal-dual active set algorithm                                   %
%                                                                         %
% INPUTS:                                                                 %
%     X   ---  matrix, normalized                                         %
%     y   ---  data vector                                                %
%    opts --- structure contains                                          %
%       lam  ---  regualrization paramter for sparsity                 %
%       tau  ---  concavity parameter                                  %
%        mu  ---  stopping parameter                                   %
%        MaxIt  ---  maximum number of iterations (default: 5)            %
%         x0    ---  initial guess for x (default: 0)                     %
%         d0    ---  initial guess for d (default: Xty)                   %
%       weigth  ---  weight* x + (1-weight)*d (default:0.5)               %
% OUTPUTS:                                                                %
%      x  ---  solution                                                   %
%      d  ---  dual                                                       %
%      G  ---  X(:,A)'*X(:,A) + alpha I                                   %
%      A  ---  active set                                                 %
%      s  ---  size of active set                                         %
%      L  ---  loss - norm(y)^2                                           %
%      It ---  # of iter for the PDAS stop                                %
%-------------------------------------------------------------------------%
n= opts.n;
p= opts.p;
implicit = isa(X,'function_handle');
lam = opts.lam;
tau = opts.tau;
MaxIt = opts.MaxIt;
weight = opts.weight;
mu = opts.mu;
Xty= opts.Xty;
Ao = opts.A;
Go = opts.G;
G = [];
so = length(Ao);
x = x0;
d = d0;
T = (2-tau)*(2*(1-tau))^((tau-1)/(2-tau))*lam^(1/(2-tau));
if implicit == 0  % initializing ...
 pd  = weight*x + (1- weight)*d;
A = find(abs(pd)>T);  %active set
s = length(A);
pdt = abs(pd)>T;  
L = 0;
It = 0;
 while  It < MaxIt && s < mu 
    It = It + 1;
    if weight == 0.5
       pd = pd/weight;
    end
     dA = zeros(s,1);  
     xA = x(A);
     Tstar = (2*lam*(1-tau))^(1/(2-tau));
     ii = find(abs(xA)>=Tstar);
     dA(ii) = lam*tau*(abs(xA(ii))).^tau./(xA(ii)); 
     XtyA = Xty(A);
    rhs =  XtyA - n*dA;
    if  s == so 
            if s == 0 
                    G = Go;
            else
                if Ao == A
                    G = Go;
                else
                    Xa = X(:,A);
                    G = Xa'*Xa ; 
                end
            end
    else
            Xa = X(:,A);
            G = Xa'*Xa ; 
    end
    x = zeros(p,1);
    xA = G\rhs;
    x(A) = xA;
    L =  dot(xA,XtyA+dA);
    sx = sparse(x);
    d = Xt*(y - X*sx)/n;
   pd  = weight*x + (1- weight)*d;
    Ao = A;
    Go = G;
    so = s;
    A = find(abs(pd)>T);
    s = length(A);
    tpd = pdt;
    pdt = abs(pd)>T;
    if length(A)> mu
        % too many nonzeros
        break;
    end
    if sum(pdt == tpd) == p
        % active set coincides
        break;
    end 
 end
  if It>0
    s = so; 
    A = Ao;
  end
else
 pd  = weight*x + (1- weight)*d;
A = find(abs(pd)>T); 
s = length(A);
pdt = abs(pd)>T;  
L = 0;
It = 0;
while  It < MaxIt && s < mu
    It = It + 1;
    if weight == 0.5
       pd = pd/weight;
    end
     dA = zeros(s,1); 
     xA = x(A);
     Tstar = (2*lam*(1-tau))^(1/(2-tau));
     ii = find(abs(xA)>=Tstar);
     dA(ii) = lam*tau*(abs(xA(ii))).^tau./(xA(ii));    
    XtyA = Xty(A); 
    rhs =  XtyA - n*dA;
     if  s == so 
            if s == 0 
                    G = Go;
            else
                if Ao == A
                    G = Go;
                else
                   option.subset = A;
                   option.p = p;
                   option.init = xA;
                   option.itcg = min(s,2); 
                   xA= Subcg(X,Xt, rhs,option);
                end
            end
    else
           option.subset = A;
           option.p = p;
           option.init = xA;
           option.itcg = min(s,2); 
           xA= Subcg(X,Xt, rhs,option);
     end      
    L =  dot(xA,XtyA+dA);
    x = zeros(p,1);
    size(x);
    x(A)=xA;    
     zz= Xt(X(x));
   d = (Xty - zz)/n;
   pd  = weight*x + (1- weight)*d;
    Ao = A;
    Go = G;
    so = s;
    A = find(abs(pd)>T);
    s = length(A);
    tpd = pdt;
    pdt = abs(pd)>T;
    if length(A)> mu
        % too many nonzeros
        break;
    end
    if sum(pdt == tpd) == p
        % active set coincides
        break;
    end 
end
if It>0
    s = so; 
    A = Ao;
end

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
%             .p         -  number of whole conlums                      %
%             .init      -  initial value                                %
%             .itcg      -  number of iterations                         %
% Output:                                                                %
%            x -  soultion                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = option.p;
maxiter = option.itcg;
subset = option.subset;
Ssubt = @(z) upsam(z,subset,p);
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


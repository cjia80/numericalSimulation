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
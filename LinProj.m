function x = LinProj(x,l1,l2,stepsize,maxoutiter,maxiniter,tolin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min_x - <Psi'y,x>/m, s.t. \|x\|_1 <= sqrt{s}, \|x\|_2<=1
% by projected gradient methods, the above model is proposed by 
%   Yaniv Plan and Roman Vershynin  "Robust 1-bit compressed sensing and 
%   sparse logistic regression: A convex programming approach", IEEE
%   Transactions on Information Theory, 59 (2013), pp. 482-494.
% by Yuling Jiao 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
grad = - x;
for k = 1:maxoutiter
    xt = x - stepsize*grad;
    x = adprojetion(xt,l1,l2,maxiniter,tolin);
end
x = x/norm(x);
end
function xlp = adprojetion(x,l1,l2,maxiniter,tolin)
%----------------------------------------------------------------
% find the solution: min_v\|v-x\|, s.t. \|x\|_1<=l1,\|x\|_2<=l2
%----------------------------------------------------------------
for k = 1:maxiniter
    tx =  projetionball(x, l1,1); % on l1 ball with radias l1
    x =   projetionball(tx, l2,2); % on l1 ball with radias l2
    if norm(x,1)<= l1+tolin && norm(x) <= l2+tolin
        break
    end
end
xlp = x;
end
function w = projetionball(v, r,mode)
%   solve the following constrained minimization problem:
%
%    min   ||w - v||_2
%    s.t.  ||w||_1 <= b, or ||w||_2 <= r
%------------------------------------------------------------------
%   Yuling Jiao 9,2016 yulingjiaomath@whu.edu.cn
%------------------------------------------------------------------------
if (r < 0)
  error('Radius of L1 ball is negative: %2.3f\n', b);
end
if mode == 1 % on l1 ball
    if (norm(v, 1) < r)
      w = v;
      return;
    end
    u = sort(abs(v),'descend');
    sv = cumsum(u);
    rho = find(u > (sv - r) ./ (1:length(u))', 1, 'last');
    theta = max(0, (sv(rho) - r) / rho);
    w = sign(v) .* max(abs(v) - theta, 0);
else
    w =  r *v / max(norm(v,'fro'),r) ;
end
end
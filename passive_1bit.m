function [xb,lambda]=passive_1bit(y,Phi,x0,K,Phit)  
 %%  
% passive_1bit is the function using the L1 penalty with the linear loss to decode the 1-bit measurements
% the code is based on \cite{yan2018nonconvex analytical solution} and
% \cite{zhang2014efficient}
% Modified by Cui Jia to work when A is a function handle.
% Modified by Cui Jia to obtain solutions with more accurate support set.
%%

 implicit = isa(Phi,'function_handle');
if implicit == 0
     if ~exist('Phit','var')
        Phit = Phi';
    end
   [m,n] = size(Phi);
   PtY = Phit*y/m;
else
  if ~exist('Phit','var')
        disp('error, funtion handle Phit is not defined')
    else
        m = length(y);
        PtY = Phit(y)/m;
        n = length(PtY);
   end
end

 least_err=1000000;
 lambda_set= (0.2:0.1: 4.1) * sqrt(log(n)/m);
 for i=1:length(lambda_set)
     lambda_temp=lambda_set(i);
     x_temp=passive(lambda_temp,PtY,n);
     temp_err = norm(x_temp - x0);
      if temp_err < least_err
            best_lambda =lambda_temp;
            least_err = temp_err;
      end
 end
lambda=best_lambda;
x=passive(lambda,PtY,n);

suppb=find(x);
if length(suppb)<K
       xb=x;
   else
       [a,b]=sort(abs(x),'descend');
       bb=b(1:K);
       xa=x(bb);
       xb=zeros(n,1);
       xb(bb)=xa;
end   
end

function x=passive(lambda_temp,PtY,n)
  v= PtY;
 lambda =lambda_temp;
 
 shrinkage_temp = lambda*ones(n,1);
 t = max(abs(v) - shrinkage_temp, 0).*sign(v);
  t_sq = t'*t;
   if t_sq > 0
     mu = sqrt(t_sq);
     x = 1/mu*t;
   else
       mu = 0;
       x = 0*t;
    end
end      
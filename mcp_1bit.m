function   [xb,par] = mcp_1bit(y, Phi, x0,K,Phit)
%%  
% mcp_1bit is the function using the mcp penalty with the linear loss to decode the 1-bit measurements
% the code is based on \cite{yan2018nonconvex analytical solution}
% Modified by Cui Jia to work when A is a function handle.
% Modified by Cui Jia to obtain solutions with more accurate support set
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

least_err = 1000000;
 mu_Set = (0.01 :0.03: 1.1) * sqrt(log(n)/m);
 b_2_Set = [0.5, 1, 3, 5, 10];
 for i = 1 : length(mu_Set)
     mu_temp = mu_Set(i);
    for j = 1 : length(b_2_Set)
        b_2_temp = b_2_Set(j); 
        x_temp = mcp_test([mu_temp,b_2_temp],PtY,m,n);
        temp_err = norm(x_temp - x0);
         if temp_err < least_err
            best_mu = mu_temp;
            best_b_2 = b_2_temp;          
            least_err = temp_err;
        end
    end
end
par = [best_mu, best_b_2];
 x = mcp_test(par,PtY,m,n);
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
function x = mcp_test(par,PtY,m,n)

       v = PtY;
       flag = 0;
        lambda = par(1);
        b = par(2);
        mu = 1/b;
        
        v_s = sort(abs(v));
        
        K = find(v_s > lambda);
        if isempty(K)
            K = m ; %delete +1
        else
            K = K(1) - 1;
        end
        
        d2 = sum(v_s(K+1 : n).^2);
        d_max = b^2*d2;
        
        if d_max > 1
%             tmp = d_max
            ii = K + 1;
            d1 = 0;
            
            while d_max > 1 && ii <= n

                mu = v_s(ii)/(b*lambda);
                d_max = d1/(mu-1/b)^2 + d2/mu^2;             
                d1 = d1 + (v_s(ii) - lambda)^2;
                d2 = d2 - v_s(ii)^2;
                d_max2 = d1/(mu-1/b)^2 + d2/mu^2;
                if d_max2 <= 1 && d_max > 1
                    flag = 1;
                    break
                end
               ii = ii + 1;
            end
%             if ii > n
%                 flag = 1;
%             end
            
%             if ii <= n || d_max <= 1
%                 d1 = d1 - (v_s(ii-1) - lambda)^2;
%                 d2 = d2 + v_s(ii-1)^2;
%             end
%             mu = sqrt(d2);
            if flag ~= 1
                mu_last = 0;
                iter_inner = 0;
                while abs(mu - mu_last) > 1e-7 && iter_inner < 200
                    mu_last = mu;
                    mu = sqrt(d1 * mu^2/(mu - 1/b)^2 +d2); 
                    iter_inner = iter_inner + 1;
                end
            end
%             d1/(mu-1/b)^2 + d2/mu^2
        else
            ii = K;
            while d_max < 1
                mu = v_s(ii)^2/(b * lambda^2);
                d_max = d2/mu^2;
                d2 = d2 + v_s(ii)^2;
                d_max2 = d2/mu^2;
                if d_max2 > 1
                    flag = 1;
                    break
                end
                ii = ii - 1;
            end
            if flag ~= 1
                d2 = d2 - v_s(ii+1)^2;
                mu = sqrt(d2);
            end
        end
        
        if mu > 1/b
            indx_1 = find(abs(v) <= lambda);
            indx_3 = find(abs(v) > b*lambda*mu);
            indx_2 = setdiff(1:n, [indx_1; indx_3]);
            x = zeros(n, 1);
            if ~isempty(indx_2)
                x(indx_2) = (abs(v(indx_2)) - lambda)/(mu - 1/b).*sign(v(indx_2));
            end
            if ~isempty(indx_3)
                x(indx_3) = v(indx_3)/mu;
            end
%             if norm(x) < 0.99
%                 mu
%                 keyboard
%             end
        else
            indx = find(v.^2 > b*lambda^2*mu);
            x = zeros(n, 1);
            x(indx) = v(indx)/mu;    
%             if norm(x) <0.99
%                 keyboard
%             end
        end
end
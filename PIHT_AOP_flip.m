function [x, iter, Loc] = PIHT_AOP_flip(b, A, x, K, L, miter, maxiter, alpha, tau,At) 
%%
% PIHT_AOP_flip is the function used to recover signals from 1-bit measurements in 
% 1-bit compressive sensing framework. It is based on the pinball loss minimization. 
% tau gives the slope of the pinball loss for the negative side. For other part
% and the parameters, please refer to BIHT_AOP_flip  
% Modeified by Yuling Jiao to work when A is a function handle
%%

implicit = isa(A,'function_handle');
if implicit == 0
    if ~exist('At','var')
        At = A';
    end
else
    if ~exist('At','var')
        disp('error, funtion handle Phit is not defined')
    end
end
    htol    = L;
    hd      = Inf;
    hd1     = Inf;
    iter    = 0;
    
    % initial "correct" location
    Loc     = 1:length(b);
    % initial "corrected" measurements
    bn      = b;
    while(htol < hd) && (iter < maxiter)  
        % x-update with PITH
        tau = tau  * 0.95;         % decrease tau 
        for inniter = 1:miter
            if implicit == 0
                 f_value = bn.*(A*x);                 % analog value
            else
                 f_value = bn.*(A(x));                 % analog value
            end
            gradient = -2*bn;                    % gradient for f_valuev < bias
            indx_p = find(f_value > 1);
            gradient(indx_p) = 2*tau*bn(indx_p); % gradient for strongly correctly classified points
            
            indx_0 = find(f_value == 1);
            gradient(indx_0) = -bn(indx_0);      % subgradient for f_value = bias
            if implicit == 0
                g = At*gradient;
            else
                g = At(gradient);
            end
            x = x - alpha * g;
            
            [~, index]  = sort(abs(x), 'descend'); % Best K-term (threshold)
            x(index(K+1:end)) = 0;
        end
         if implicit == 0
                y_t  = A * x;
         else
                y_t  = A(x);
         end
        
        hd   = nnz(b - sign(y_t));

       if hd < hd1 % only update Loc and bn when fewer sign flips are found
            hd1     = hd;
            % find the largest L elements of phi(b,Ax)
            [~, index] = sort(abs(y_t) .* max(-sign(b.*y_t), 0), 'descend');
            % update the "correct" location
            Loc     = 1:length(b);
            Loc(index(1:L)) = [];
            % correct the measurements by flipping the sign of b at these
            % locations
            bn      = b;
            bn(index(1:L)) = -bn(index(1:L));
            alpha   = 1;
        end
        iter = iter + 1;
    end
end


function x = BIHT_1(y, Phi, K,Phit)
% Modeified by Yuling Jiao to work when Phi is a function handle
implicit = isa(Phi,'function_handle');
if implicit == 0
    N  = size(Phi,2);
    A = @(in) sign(Phi*in);
    if ~exist('Phit','var')
        Phit = Phi';
    end
else
    if ~exist('Phit','var')
        disp('error, funtion handle Phit is not defined')
    end
    N = length(Phit(y));
    A = @(in) sign(Phi(in));
end

if K >=100
   maxiter = 200;
else
   maxiter = 100;
end
htol = 0;

x = zeros(N,1);
hd = Inf;

ii=0;
while(htol < hd)&&(ii < maxiter)
	% Get gradient
    if implicit == 0
	   g = Phit*(A(x) - y);
    else
       g = Phit(A(x) - y);
    end
	% Step
	a = x - g;

	% Best K-term (threshold)
	[trash, aidx] = sort(abs(a), 'descend');
	a(aidx(K+1:end)) = 0;

	% Update x
	x = a;

	% Measure hammning distance to original 1bit measurements
	hd = nnz(y - A(x));
	ii = ii+1;
end
% Now project to sphere
x = x/norm(x);












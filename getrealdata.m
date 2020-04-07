
function [X,Xt,y,xe,A,out] = getrealdata(imgsize,nsam,sigma,seednum)  
    rand('state',seednum)
    randn('state',seednum)
    N = imgsize;
    h = 1/N;
    mesh = (h/2:h:1-h/2)';
    out.mesh = mesh;
    Img  = -8*(.25*(.10<mesh&mesh<.102))+3*(.75*(.29<mesh&mesh<.292))... 
     + 10*(.25*(.57<mesh&mesh<.572))+ ((.75<mesh&mesh<.753).*-2);
    Img = sign(Img);
    Imgbar = mean(Img(:));
    out.bar = Imgbar;
    Img = Img - Imgbar;
   %  Sampling operator
    R = randn(nsam,N);
    Phi = @(z) R*z;     
    Phit = @(z) R'*z;   
% Effective sensing operator
   dW_L = 1;                    % levels of wavelet transform
   wav = 'db1';                 % type of wavelet 
   [lo,hi,rlo,rhi] = wfilters(wav);
   W = @(z) WaveDecRec(z,dW_L,lo,hi,rlo,rhi,1,1,0);  % Rec
   Wt = @(z) WaveDecRec(z,dW_L,lo,hi,rlo,rhi,1,1,1); % Dec
   xe =  Wt(Img);
   id = find(abs(xe)<= 1e-1);
   xe(id) = zeros(length(id),1);
   xe = xe/norm(xe);
   A = find(xe);
   Imr = (W(xe)+ Imgbar)/norm(xe); 
   figure(1)
   plot(mesh,Imr(:),'k-','LineWidth',1)
   axis off
   box off
   title('True signal');
   noise = sigma*randn(nsam,1);
   X = @(z) Phi(Wt(z));
   Xt = @(z) W(Phit(z));
   y = X(xe) + noise;
   y = sign(y);
   out.I = Imr;
   out.W = W;
end

function [Y,strc] = WaveDecRec(X,level,lo,hi,rlo,rhi,emd,ds,isdec)
persistent s;
if isdec == 1
   [Y,s] = mwavedec(X,level,lo,hi,emd,ds);
    strc = s;
    Y = Y';
else
    Y  = mwaverec(X,s,rlo,rhi,emd,0);
    strc = s;
end
end

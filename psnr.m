function y = psnr(im1,im2,scale)
if nargin< 3
    scale = 1;
end
[m,n] = size(im1);
x1 = double(im1(:));
x2 = double(im2(:));
sq_err = norm(x1-x2);
sq_err = (sq_err*sq_err)/(m*n);
if sq_err > 0
    y = 10*log10(scale^2/sq_err);
else
    disp('infinite psnr');
end
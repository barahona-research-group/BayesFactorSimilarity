function [yym,hyp2] = gpRegression(y,x,hyp2)
% 1-D regression 
% y matrix size N*n, N samples
% x vector length n

[N,n] = size(y);

xx = linspace(min(x), max(x), 100)';
% xx = log(linspace(0,10,1000)'+1);

x = repmat(x(:)',N,1);

meanfunc = @meanConst; %xmean: the model from x to f f= gp(xmean,K)
hyp.mean = mean(y(:));

likfunc = @likGauss;
sn = exp(0); 
hyp.lik = log(sn); %the model from f to Y, noise level

covfunc = @covSEiso; % covariance function
% ell = exp(-0.3523);
ell = 2*(max(x(:))- min(x(:)))/n;
% sf = exp(-0.4399);
sf = 1;
hyp.cov = [log(ell); log(sf)];  %K: the model from x to f f= gp(xmean,K)

if nargin == 2
    hyp2 = hyp;
    hyp2 = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x(:), y(:));
end
[yym yys fm fs] = gp(hyp2, @infExact, meanfunc, covfunc, likfunc, x(:), y(:), xx); %prediction
vis = 1;
if vis ==1
    figure
    f = [yym+2*sqrt(yys); flipdim(yym-2*sqrt(yys),1)];
    hold on;
    fill([xx; flipdim(xx,1)], f, [7 7 7]/8)
    plot(x, y, 'b*');plot(xx, yym,'r');hold off
end

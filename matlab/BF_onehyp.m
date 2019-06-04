function [bfmatrix] = BF_onehyp(data,x,hyp_in)

%using GPML toolbox http://www.gaussianprocess.org/gpml/code/matlab/doc/

%Gaussian Process Regression likelihood ratio matrix
%
% [bfm,likv]  = BF2(x,data)
% x: time points as a column vector
% data: data matrix. each row is one sample time series.
%
% bfmatrix: Similarity Matrix.  bfmatrix(i,j) = log( p(yi,yj) / (p(yi)p(yj)) )


assert(size(data,1)==length(x));

n = size(data,1);
bfmatrix = zeros(n,n);

meanfunc = @meanConst; %mean function
hyp.mean = mean(data(:));

likfunc = @likGauss; %likelihood function
covfunc = @covSEiso; %covariance function
sf = 1;
ell = (max(x(:))- min(x(:)))/(length(unique(x)));
hyp.cov = [log(ell); log(sf)];% initial values log(l), log(sf)
sn = 0.1
hyp.lik = log(sn); %noise level

if nargin == 3  % the hyper parameters are given
    hyp = hyp_in;
else
    % learn the hyper parameters
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x, data');
end

x12 = [x;x];
% likelihood of two data points
for i = 1:n-1
    y1 = data(i,:)'; 
    for j = i+1:n  %parfor
        y2 = data(j,:)';
        y12 = [y1;y2];
        nlik_i = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x, [y1,y2]); %negative log likelihood of y1,y2 
        nlik12 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x12, y12); %negative log likelihood of y1;y2 
        bfmatrix(i,j) = nlik_i - nlik12;   
        
        %[~,~,~,~,lp_i] = gp(hyp, @infExact, [], covfunc, likfunc, x, [y1,y2],x,[y1,y2]); %log predictive prob of y1,y2 
        %[~,~,~,~,lp_12] = gp(hyp, @infExact, [], covfunc, likfunc, x12, y12, x12,  y12); %log predictive prob of [y1;y2] 
        %bf_post(i,j) = sum(lp_12) - sum(lp_i(:)); 

    end
    %fprintf('i= %4i;\n', i);
end

bfmatrix =  bfmatrix + bfmatrix';
bfmatrix(1:n+1:end) = 0;

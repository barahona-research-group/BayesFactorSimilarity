function [bfmatrix,bf_post ]= BF_twohyp(data,x,hyp)
% Compute the GP based similarity when the time courses are synchronous 
% For a pair of time courses, optimise the hyperparameters for each time 
% course seperately.
% 
% using GPML toolbox http://www.gaussianprocess.org/gpml/code/matlab/doc/

%Gaussian Process Regression likelihood ratio matrix
%
% [bfm,likv]  = BF_twohyp(x,data)
% x: time points as a column vector
% data: data matrix. each row is one sample time series.
%
% bfmatrix: Similarity Matrix.  bfm(i,j) = log( p(yi,yj) / (p(yi)p(yj)) )

if nargin == 3  % TODO: do not learn the hyper parameters
end

n = size(data,1);
bfmatrix = zeros(n,n);
bf_post  = zeros(n,n);

meanfunc = @meanConst; %xmean: the model from x to f f= gp(xmean,K)
hyp.mean = mean(data(:));
likfunc = @likGauss;
covfunc = @covSEiso;
hyp.cov = [0; 0];% initial values log(l), log(sf)
hyp.lik = log(4);

hyp2 = hyp;
hyp2.lik = log(4);
x12 = [x;x];

for i = 1:n-1
    y1 = data(i,:)'; 
    parfor j = i+1:n  %parfor
        y2 = data(j,:)';
        y12 = [y1;y2];
        
        hyp_i = minimize(hyp, @gp, -100, @infExact, [], covfunc, likfunc, x, [y1,y2]);
        
        %hyp_i.lik = log(0.000001);
        nlik_i = gp(hyp_i, @infExact, [], covfunc, likfunc, x, [y1,y2]); %negative log likelihood of y1,y2 

        hyp12 = minimize(hyp2, @gp, -100, @infExact, [], covfunc, likfunc, x12, y12);
        
        %hyp12.lik = log(0.000001);
        nlik12 = gp(hyp12, @infExact, [], covfunc, likfunc, x12, y12); %negative log likelihood of y1;y2 
        bfmatrix(i,j) = nlik_i - nlik12;
        
        [~,~,~,~,lp_i] = gp(hyp_i, @infExact, [], covfunc, likfunc, x, [y1,y2],x,[y1,y2]); %log predictive prob of y1,y2 
        [~,~,~,~,lp_12] = gp(hyp12, @infExact, [], covfunc, likfunc, x12, y12, x12,  y12); %log predictive prob of [y1;y2] 
        bf_post(i,j) = sum(lp_12) - sum(lp_i(:)); 
    end
    % fprintf('i= %4i;\n', i);
end

bfmatrix =  bfmatrix + bfmatrix';
bfmatrix(1:n+1:end) = 0;

bf_post = bf_post + bf_post';
bf_post(1:n+1:end) = 0;

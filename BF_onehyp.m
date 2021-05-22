function [bfmatrix,bf_bregman] = BF_onehyp(data,x,varargin)
% Compute the GP based similarity when the time courses are synchronous
%using GPML toolbox http://www.gaussianprocess.org/gpml/code/matlab/doc/

%Gaussian Process Regression likelihood ratio matrix
%
% [bfm,likv]  = BF_onehyp(x,data,nRep)
% x: time points as a column vector
% data: data matrix. each row is one sample time series. 
%       in the form of [y(1,t1),...y(1,tT),y(2,t1),...,y(2,tT),...y(nRep,t1),...y(nRep,tT)]
% nRep: number of biological replicates.
%
% bfmatrix: Similarity Matrix.  bfm(i,j) = log( p(yi,yj) / (p(yi)p(yj)) )

if nargin < 4   
    nRep = 1; 
else
    nRep = varargin{2}; %biological replicates.
end

n = size(data,1);
nTime = length(x);
bfmatrix = zeros(n,n);
bf_bregman = zeros(n,n);

meanfunc = @meanConst; %xmean: the model from x to f f= gp(xmean,K)
hyp.mean = mean(data(:));
likfunc = @likGauss;
covfunc = @covSEiso;
hyp.cov = [log(1); log(1)];% initial values log(l), log(sf)
hyp.lik = log(2); %noise level

% learn the hyper parameters
x1 = repmat(x,nRep,1);
x12 = [x1;x1];

if nargin == 3  % do not learn the hyper parameters
    hyp = varargin{1};
    disp('hyper-parameter provided.')
    disp(hyp)
else
    hyp = minimize(hyp, @gp, -100, @infExact, meanfunc, covfunc, likfunc, x1, data');
    disp('hyper-parameter learned.')
    disp(hyp)
end
% likelihood of two data points
% sn2 = exp(2*hyp.lik); W = ones(nTime,1)/sn2;            % noise variance of likGauss
% K = apx(hyp,{covfunc},x1,[]);                        % set up covariance approximation
% [ldB2,solveKiW,dW,dhyp,post.L] = K.fun(W); % obtain functionality depending on W
% 
% alpha = solveKiW(data' - hyp.mean);
% V = K.P(alpha);  
K = feval(covfunc, hyp.cov, x1);Ky = eye(nTime) * exp(2*hyp.lik) + K;
V1 = Ky\(data'-hyp.mean);V = V1;

for i = 1:n
    y1 = data(i,:)'; 
    for j = i:n  %parfor
        y2 = data(j,:)';
        y12 = [y1;y2];
        nlik_i = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x1, [y1,y2]); %negative log likelihood of y1,y2 , independnet
        nlik12 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x12, y12); %negative log likelihood of y1;y2 
        bfmatrix(i,j) = nlik_i - nlik12;
        bf_bregman(i,j) = V(:,i)'*K*V(:,i) + V(:,j)'*K*V(:,j) - 2* V(:,i)'*K*V(:,j);
%        [~,~,~,~,lp_i] = gp(hyp, @infExact, [], covfunc, likfunc, x1, [y1,y2],x,[y1,y2]); %log predictive prob of y1,y2 
%        [~,~,~,~,lp_12] = gp(hyp, @infExact, [], covfunc, likfunc, x12, y12, x12,  y12); %log predictive prob of [y1;y2] 
%        bf_post(i,j) = sum(lp_12) - sum(lp_i(:)); 
    end
    fprintf('i= %4i;\n', i);
end


bfmatrix =  bfmatrix + bfmatrix';
%bfmatrix(1:n+1:end) = 0;

bf_bregman =  bf_bregman + bf_bregman';
%bf_post = bf_post + bf_post';
%bf_post(1:n+1:end) = 0;

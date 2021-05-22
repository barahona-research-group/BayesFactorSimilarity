function [bfmatrix,bf_bregman] = BF_async(data,x,mask,varargin)
% Compute the GP based similarity when the time courses are asynchronous
% Use a mask matrix to represent missing data.
% using GPML toolbox http://www.gaussianprocess.org/gpml/code/matlab/doc/

%Gaussian Process Regression likelihood ratio matrix
% 
% [bfm,likv]  = BF_async(x,data,nRep)
% x: time points as a column vector
% data: data matrix. each row is one sample time series. n X nTime
% mask: mask matrix of missing data of 0 and 1

% bfmatrix: Similarity Matrix.  bfm(i,j) = log( p(yi,yj) / (p(yi)p(yj)) )

nRep = 1; 

if nargin == 4  % TODO: do not learn the hyper parameters
end

n = size(data,1);
nTime = size(x,1);
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

if nargin == 4  % do not learn the hyper parameters
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
    id_i = find(mask(i,:));
    y1 = data(i,id_i)';
    x_i = x1(id_i);
    v_i = Ky(id_i,id_i)\y1;
    
    nlik_i = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x_i, y1); %negative log likelihood of y1,y2 
    for j = i:n  %parfor      
        id_j = find(mask(j,:));
        y2 = data(j,id_j)';
        x_j = x1(id_j);
        y12 = [y1;y2];
        x_ij = [x_i;x_j];
        nlik_j = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x_j, y2); %negative log likelihood of y1,y2 
        nlik12 = gp(hyp, @infExact, meanfunc, covfunc, likfunc, x_ij, y12); %negative log likelihood of y1;y2 
        bfmatrix(i,j) = nlik_i + nlik_j - nlik12;
        
        v_j = Ky(id_j,id_j)\y2;

        bf_bregman(i,j) = v_i'*K(id_i,id_i)*v_i ...
            + v_j'*K(id_j,id_j)*v_j - 2* v_i'*K(id_i,id_j)*v_j;

    end
    fprintf('i= %4i;\n', i);
end


bfmatrix =  bfmatrix + bfmatrix';
%bfmatrix(1:n+1:end) = 0;

bf_bregman =  bf_bregman + bf_bregman';
%bf_post = bf_post + bf_post';
%bf_post(1:n+1:end) = 0;

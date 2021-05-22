function [c,vi] = partition_combinatorial(W,Nc)
% W the adjacency matrix
% Nc number of clusters
% not for large sparse 

% output Discrete, the indicator matrix H.
% Zijing Liu, Dec 2015


%TODO
%check input 

%check symmetry


L = diag(sum(W)) - W; %combinatorial Laplacian

n = size(W,1);
   dataNcut.offset = 5e-1;
    dataNcut.verbose = 0;
    dataNcut.maxiterations = 300;
    dataNcut.eigsErrorTolerance = 1e-8;
    dataNcut.valeurMin=1e-6;
    options.issym = 1;
     
if dataNcut.verbose
    options.disp = 3; 
else
    options.disp = 0; 
end
options.maxit = dataNcut.maxiterations;
options.tol = dataNcut.eigsErrorTolerance;

options.v0 = ones(size(L,1),1);
options.p = max(35,2*Nc); %voir
options.p = min(options.p,n);

L = sparsifyc(L,dataNcut.valeurMin);
% [v,d,convergence] = eigs(@mex_w_times_x_symmetric,size(L,1),Nc,'sa',options,tril(L)); 

[v,d] = eigs(L,Nc,'sa',options); %computation load:eigs


d = diag(d);

[d,is] = sort(d);
v = v(:,is);
Eigenvectors = v;

parfor i = 1:100
    [Discrete,cEigenvectors] =discretisation(Eigenvectors);
    Discrete = full(Discrete);
    C(:,i) = transformHMatrixtoPartitionVector(Discrete);
end
[Stb,Nc,c] = stability_maxcb(W,C,1); %computation load: 100 expm 
vi = varinfo(C',true);

% Discrete = full(Discrete);
function [Cm,Ch,C_ncut,C_njw] = cluster_dist_matrix(sm,nc,data,dm,k)
% Clustering using the similarity matrix
% input -- sm: the similarity matrix
%          nc: number of communities to find
%          data: the data matrix
%          dm: distance matrix 

n = size(sm,1);

if nargin == 3
    dm = 1 + max(sm(:)) - sm;
    dm(1:n+1:end) = 0;
    k = 7;
end

if nargin == 4
    k = 7;
end
Cm = []
Z = linkage(squareform(dm),'average');
Ch = cluster(Z,'maxclust',nc);

neighbor_num = k;
G = constructNetworkStructure(data',dm,'knn',neighbor_num);
A = double(G);
[c_ncut,x] = ncutW(A,nc);
c_ncut = transformHMatrixtoPartitionVector(c_ncut);
[c1,x] = gcut(A,nc);
c_njw = c_ncut;
    for i = 1:length(c1)
        c_njw(c1{i}) = i;
    end
C_ncut = c_ncut;
C_njw = c_njw;

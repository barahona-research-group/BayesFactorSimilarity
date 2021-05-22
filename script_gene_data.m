addpath(genpath('./lib/gpml-matlab-v4.2-2018-06-11'))
addpath(genpath('./lib/Ncut_9'))
addpath(genpath('./lib/ZPclustering'))

load('gene_data.mat')

%% compute distances
[bfm,d_breg] = BF_onehyp(data_s,t);

sm = bfm;
n = length(sm);
dist_bf = 1 + max(sm(:)) - sm;
dist_bf(1:n+1:end) = 0;

dist_dtw = zeros(n,n);
for i=1:n
    for j=1:n
        dist_dtw(i,j) = dtw(data_s(i,:),data_s(j,:));
    end
end
%% hierarchical clustering
n_candi = 3:17;

lkge = 'complete';

Z = linkage(squareform(dist_bf),lkge);
ch = cluster(Z,'maxclust',n_candi);

dist_eu = squareform(pdist(data_s));
Z_eu = linkage(squareform(dist_eu),lkge);
ch_eu = cluster(Z_eu,'maxclust',n_candi);

Z_dtw = linkage(squareform(dist_dtw),lkge);
ch_dtw = cluster(Z_dtw,'maxclust',n_candi);

writeout_partition([ch,ch_eu,ch_dtw],'bf_complete_dtw.csv');

%% spectral clustering
neighbor_num = 10;

n_candi = 3:17;
D = dist_bf;
for k = 1:length(n_candi)
    n_c = n_candi(k);
    G = constructNetworkStructure(data',D,'knn',neighbor_num);
    A = double(G);
    [c_ncut,x] = ncutW(A,n_c);
    c_ncut = transformHMatrixtoPartitionVector(c_ncut);
    [c1,x] = gcut(A,n_c);
    c_njw = c_ncut;
    for i = 1:length(c1)
        c_njw(c1{i}) = i;
    end
    C_ncut(:,k) = c_ncut;
    C_njw(:,k) = c_njw;
end

D = squareform(pdist(data_s,'euclidean'));
for k = 1:length(n_candi)
    n_c = n_candi(k);
    G = constructNetworkStructure(data',D,'knn',neighbor_num);
    A = double(G);
    [c_ncut,x] = ncutW(A,n_c);
    c_ncut = transformHMatrixtoPartitionVector(c_ncut);
    [c1,x] = gcut(A,n_c);
    c_njw = c_ncut;
    for i = 1:length(c1)
        c_njw(c1{i}) = i;
    end
    C_ncut_eu(:,k) = c_ncut;
    C_njw_eu(:,k) = c_njw;
end


D = dist_dtw;
for k = 1:length(n_candi)
    n_c = n_candi(k);
    G = constructNetworkStructure(data',D,'knn',neighbor_num);
    A = double(G);
    [c_ncut,x] = ncutW(A,n_c);
    c_ncut = transformHMatrixtoPartitionVector(c_ncut);
    [c1,x] = gcut(A,n_c);
    c_njw = c_ncut;
    for i = 1:length(c1)
        c_njw(c1{i}) = i;
    end
    C_ncut_dtw(:,k) = c_ncut;
    C_njw_dtw(:,k) = c_njw;
end

D = d_breg;
for k = 1:length(n_candi)
    n_c = n_candi(k);
    G = constructNetworkStructure(data',D,'knn',neighbor_num);
    A = double(G);
    [c_ncut,x] = ncutW(A,n_c);
    c_ncut = transformHMatrixtoPartitionVector(c_ncut);
    [c1,x] = gcut(A,n_c);
    c_njw = c_ncut;
    for i = 1:length(c1)
        c_njw(c1{i}) = i;
    end
    C_ncut_breg(:,k) = c_ncut;
    C_njw_breg(:,k) = c_njw;
end

% save the partitions
writeout_partition([C_ncut_dtw,C_njw_dtw],'sc_dtw.csv');
writeout_partition([C_ncut,C_njw,C_ncut_eu,C_njw_eu],'sc_bf_eu.csv');
writeout_partition([C_ncut_breg,C_njw_breg],'sc_breg.csv');

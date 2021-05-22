addpath(genpath('/home/zliu2/Codes/PartitionStability'))
addpath(genpath('/home/zliu2/Codes/RBS'))
addpath(genpath('/home/zliu2/Codes/MarkovZoomingMap'))
addpath(genpath('/home/zliu2/Codes/Zijing'))
addpath(genpath('/home/zliu2/Codes/Ncut_9'))

A = W;%graph matrix

tic;
nc = 100; %max #groups

C_rcut = zeros(length(A),nc);
vi = zeros(1,nc);
for c = 2:nc
    nbCluster = c;
    [C_rcut(:,c),vi(c)] = partition_combinatorial(A,nbCluster);
end
toc;

tic;
% Time = 10.^[0:0.1:4];
Time = logspace(-3,5,600);
[Stb,Nc,C] = stability_maxcb(A,C_rcut,Time);
toc;
%stability_plot(Time,Stb,Nc,vi(Nc))

save('spectral_rcut.mat','C','Time','Stb','Nc','vi');

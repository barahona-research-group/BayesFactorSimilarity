run('./lib/gpml-matlab-v4.2-2018-06-11/startup')
addpath(genpath('./lib/InfoTheory'));
addpath(genpath('./lib/PartitionStability'));
addpath(genpath('./lib/Ncut_9'))
addpath(genpath('./lib/ZPclustering'))

% test the BF similarity for time course data.
% generate data
%%
n_t = 15; %length of the time course
sigma = 0.15; % level of noise
n = 150; % number of the time course
n_g = 3; % 3 clusters
c_g = zeros(150,1);c_g(51:100) = 1;c_g(101:end) = 2; %cluster ids

t1 = linspace(0,1,n_t)'; % regular sampling
% t1 = [linspace(0,0.25,10)';linspace(0.3,1,5)']; % irregular sampling

n_run = 10; % repeat 10 times
result_fname = ['rslt_',num2str(sigma*100),'e.mat'];
%% 
% the output NMI values

% hierarchical clustering
rh_eu = zeros(n_run,1); %Eulidean 
rh_dtw = zeros(n_run,1); %DTW 
rh_bf = zeros(n_run,1); % Proposed Bayes Factor
rh_breg = zeros(n_run,1); % Bregman

% spectral clustering
rsp_eu = zeros(n_run,1);
rsp_dtw = zeros(n_run,1);
rsp_bf = zeros(n_run,1);
rsp_breg = zeros(n_run,1);

%
parfor n_iter = 1:n_run
    disp(n_iter)
    x1 = zeros(n_t,n);
    for i = 1: n/3
        x1(:,i) = t1 + 0.25*(t1).^4 + sigma.*randn(n_t,1);
        x1(:,i+n/3) = t1 + sigma.*randn(n_t,1);
        x1(:,i+n/3*2) = t1-0.25*(t1).^4 + sigma.*randn(n_t,1);
    end
    x1 = x1-0.5;
    
    d_eu = squareform(pdist(x1'));
    [s_bf,d_breg] = BF_onehyp(x1',t1);
    
    [Cm_bf,ch_bf,c_ncut_bf] = cluster_dist_matrix(s_bf,n_g,x1);
    [Cm_eu,ch_eu,c_ncut_eu] = cluster_dist_matrix(-d_eu,n_g,x1,d_eu);
    [Cm_breg,ch_breg,c_ncut_breg] = cluster_dist_matrix(-d_breg,n_g,d_breg);
    
    %spectral clustering
    rsp_breg(n_iter) = nmi(c_g,c_ncut_breg);
    rsp_eu(n_iter) = nmi(c_g,c_ncut_eu);
    rsp_bf(n_iter) = nmi(c_g,c_ncut_bf);
    %hierachical clustering
    rh_breg(n_iter) = nmi(c_g,ch_breg);
    rh_eu(n_iter) = nmi(c_g,ch_eu);
    rh_bf(n_iter) = nmi(c_g,ch_bf);
    
end
%% dtw
for n_iter = 1:n_run
    dist_dtw = zeros(n,n);
    for i =1:n
        for j =1:n
            dist_dtw(i,j) = dtw(x1(:,i),x1(:,j));
        end
    end
        
    [Cm_dtw,ch_dtw,c_ncut_dtw] = cluster_dist_matrix(-dist_dtw,n_g,x1,dist_dtw);

    rsp_dtw(n_iter) = nmi(c_g,c_ncut_dtw);

    rh_dtw(n_iter) = nmi(c_g,ch_dtw);
     
end

rslt = [rsp_bf,rsp_eu,rsp_breg,rsp_dtw,rh_bf,rh_eu,rh_breg,rh_dtw];
% save the result
save(result_fname,'rslt')


%% plot
y = rslt;

errlow = std(y);
errhigh = std(y);
y = mean(y);
figure,hold on
bar(1:4,y(1:4),0.4,'FaceColor',[0.86,0.86,0.86])
bar(6:9,y(5:8),0.4,'FaceColor',[0.95,0.87,0.73])

er = errorbar([1:4,6:9],y,errlow,errhigh);    
er.Color = [0 0 0];                            
er.LineStyle = 'none';
ylabel('NMI')
xlim([0.5,9.5]);ylim([0,1.1])
xticks([1,2,3,4,6,7,8,9])
xticklabels({'Proposed','Euclidean','Bregman','DTW','Proposed','Euclidean','Bregman','DTW'})
text(1.8,1.05,'Spectral clustering')
text(6.6,1.05,'Hierarchical clustering')
hold off
set(gcf,'position',[100,100,600,200])
title(['noise=',num2str(sigma)])
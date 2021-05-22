function [Stb,Nc,C] = stability_maxcb(W,c,time)
% Find Markov stability of given partitions at given times
% combinatorial Laplacian
% Inputs:
% W ----- the weighted adjacency matrix N*N
% c ----- partitions (N*n) N: number of nodes; n: number of partitions
% every col is a partition vector
% time -- the Markov time vector
% Outputs:
% Stb --- the stability
% Nc ---- the number of communities
% C ----- the optimal partition for each Markov time contained in time. size N*t
%
% Zijing Liu, Dec 2015


% input check //todo
% //


[N,n] = size(c);
H = cell(1,n);
% transform to H matrix
for i = 1:n
    H{i} = full(transformPartitionVectorToHMatrix(c(:,i)));
end

% some related matrices
D = diag(sum(W));
L = D-W; % Laplacian
Ld = L./sum(W(:))*N;
P = ones(N,1)/N; %stationary distribution vector
Stb = zeros(length(time),1);
Nc = zeros(length(time),1);
C = zeros(N,length(time));
if length(time)<4  
    for i = 1:length(time)
        E = expm(-time(i).* Ld); %continuous time
        % E = (eye(n)-Ld)^time; %descrete time
        S = diag(P)*E - P*P';
        s = zeros(1,n);
        for j = 1:n
            h = H{j};
            s(j) = trace(h'*S*h); %stabilities for this time
        end
        [Stb(i),Nc(i)] = max(s); %find the partition giving the best stability
        C(:,i) = c(:,Nc(i));
        Nc(i) = max(C(:,i)) - min(C(:,i)) + 1;
%         disp(['Time',num2str(i),'done'])
    end
else %parallel
    parfor i = 1:length(time)
        E = expm(-time(i).* Ld); %continuous time
        % E = (eye(n)-Ld)^time; %descrete time
        S = diag(P)*E - P*P';
        s = zeros(1,n);
        for j = 1:n
            h = H{j};
            s(j) = trace(h'*S*h); %stabilities for this time
        end
        [Stb(i),Nc(i)] = max(s); %find the partition giving the best stability
        C(:,i) = c(:,Nc(i));
        Nc(i) = max(C(:,i)) - min(C(:,i)) + 1;
        disp(['Time',num2str(i),'done'])
    end
end

    function H = transformPartitionVectorToHMatrix(pvector)
        %Transforms a given partition vector into a valid H matrix (Partition incidence matrix)
        % assumes that the ordering is contigous and starts from either zero or
        % one.
        
        if min(pvector)==0
            pvector = pvector+1;
        end
        nr_nodes = length(pvector);
        
        % create H matrix
        H = sparse(1:nr_nodes,pvector,1);
           
    end

end
function [G, Ij, mj] = clustering_GCP(eff, groups)
%% clusetering all the segments one by one up to only one
len = length(groups(1,:));% number of segments
mj = groups(2,:)-groups(1,:)+1;% the number of data points in each segment
Ij = zeros(1,len);% the averaged intensity of each segment
G = struct([]);% the structure to store all the clustering history
for i = 1:len
    Ij(i) = mean(eff(groups(1,i):groups(2,i)));
    G(1).g(i).gg = i;
end
G = sub_clustering(Ij, mj, G);% clustering
end

%% 
% for each cycle, this sub function cluster two segments into one
function G = sub_clustering(Ij, mj, G)
n = numel(G(end).g);
if n <= 1
    return
else
    M = ones(n, n)*(-inf);% the merit matrix
    for i = 1:n-1
        for j = i+1:n
            % calculating the merit for all the possible pairs
            M(i, j) = llr_merit(Ij(i), Ij(j), mj(i), mj(j));
        end
    end
    while numel(G(end).g) > 1;
        numel(G(end).g)
        [a, b] = size(M);
        [~, q] = max(M(:));
        j = ceil(q/a);% find the column
        i = q-(j-1)*a;% find the row
        G(end+1).g = G(end).g;% initialize G(n+1)
        G(end).g(i).gg = [G(end).g(i).gg, G(end).g(j).gg];% group segments i and j together
        G(end).g(j) = [];% free another space
        Ij(i) = (Ij(i)*mj(i)+Ij(j)*mj(j))/(mj(i)+mj(j));% group Ij
        Ij(j) = [];
        mj(i) = mj(i)+mj(j);% group mj
        mj(j) = [];
        %% update the M
        M(:, i) = -inf; M(i, :) = -inf;
        M(:, j) = []; M(j, :) = [];
        for k = 1:numel(G(end).g)% pairing with the new cluster i
            if k < i
                M(k, i) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));
            elseif k > i
                M(i, k) = llr_merit(Ij(i), Ij(k), mj(i), mj(k));
            end% if
        end% for k
    end% while
end% if n<=1
end% function

%% the merit function for Gaussian model (but can be more general)
function llr = llr_merit(I1, I2, m1, m2)
I = (I1*m1 + I2*m2)/(m1+m2);
llr = (m1+m2)*I^2-m1*I1^2-m2*I2^2;
end
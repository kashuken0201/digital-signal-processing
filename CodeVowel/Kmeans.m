function Kmean = Kmeans(x, K_means)
[r, c] = size(x);
sizeCluster=ceil(r/K_means);
Kmean = zeros(K_means,c);

for i=1:r
    if mod(i,sizeCluster)==0
        Kmean(i/sizeCluster,:)=Kmean(i/sizeCluster,:)+x(i,:);
        Kmean(i/sizeCluster,:)=Kmean(i/sizeCluster,:)/sizeCluster;
    else
        Kmean(floor(i/sizeCluster)+1,:)=Kmean(floor(i/sizeCluster)+1,:)+x(i,:);
    end
    if i==r && mod(i,sizeCluster) ~= 0
        Kmean(floor(i/sizeCluster)+1,:)=Kmean(floor(i/sizeCluster)+1,:)/(mod(r,sizeCluster));
    end
end
end
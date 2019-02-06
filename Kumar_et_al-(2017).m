function [ A ] = Kumar_et_al( Imagem )
%% Title: Automated diagnosis of acute lymphocytic leukemia and acute myeloid leukemia using multi-SV
% Authors: Kumar PS, Vasuki S
% Journal of Biomedical Imaging and Bioengineering (2017) Volume 1, Issue 1

%% RGB to HSV
HSV = rgb2hsv(Imagem);
V = HSV(:,:,3);

%% K-means application

% Hue and Saturation selection. 
HS = double(HSV(:,:,1:2));
nrows = (size(HS,1));
ncols = (size(HS,2));
HS = reshape(HS, nrows*ncols,2);

% Number of clusters
Clusters = 4;

% K-means execution
[idx, center] = kmeans(HS, Clusters, 'distance', 'sqEuclidean',...
    'Replicates',5);
pixel_labels = reshape(idx, nrows, ncols);
mean_cluster_value = mean(center,2);
[tmp, idx] = sort(mean_cluster_value);
nuclei_cluster = idx(4);

% Image reconstruction
[u,v] = size(V);
for i = 1:u
    for j = 1:v
           if pixel_labels(i,j) == nuclei_cluster
              A(i,j) = 1; 
           else
              A(i,j) = 0;
           end 
    end
end
end


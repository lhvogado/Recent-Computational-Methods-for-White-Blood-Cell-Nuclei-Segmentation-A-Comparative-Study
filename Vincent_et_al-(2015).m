function [ pixel_labels ] = WBC_Vincent_et_al( I )
%% Lendo as imagens 
transforma = makecform('srgb2lab');
lab_he = applycform(I,transforma);

%% K-means
%Agrupar as classes para se adequar aos parâmetros que o K-Means requer.
ab = double(lab_he(:,:,2:3));
nrows = (size(ab,1));
ncols = (size(ab,2));
ab = reshape(ab,nrows*ncols,2);

%Três clusters
nColors = 3;

[cluster_idx, cluster_center] = kmeans(ab, nColors, 'distance','sqEuclidean',...
  'Replicates',3);%Aplico o Kmeans
pixel_labels = reshape(cluster_idx, nrows, ncols);

% imshow(pixel_labels);
%% Refazendo a imagem
mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
nuclei_cluster = idx(1);
[u,v] = size(I);
v = v/3;
for i = 1:u
    for j = 1:v
           if pixel_labels(i,j) == nuclei_cluster
              A(i,j) = 0; 
           else
               A(i,j) = 1;
           end 
    end
end
R = I(:,:,1);
G = I(:,:,2);
B = I(:,:,3);
for i = 1:u
    for j = 1:v
           if A(i,j) ~= 0
              R(i,j) = 0;
              G(i,j) = 0;
              B(i,j) = 0;
           end 
    end
end

Final = cat(3,R,G,B);
R = Final(:,:,1);
G = Final(:,:,2);
B = Final(:,:,3);
R = imadjust(R);
B = imadjust(B);
G = imadjust(G);
I = cat(3,R,G,B);

%% OTSU
level = graythresh(I);
BW = im2bw(I,level);
%% Morph
se = strel('disk', 1);
se1 = strel('disk', 5);
IM = imopen(BW, se);
IM = imclose(IM, se1);
end


%%  Gera as imagens com os algoritmos de Madhloom e de o Proposto por Mohammed

clc;
clear;
caminho = 'C:\Users\Henrique\OneDrive\Matlab\Implementações-Segmentação\LEUK_Resultados\Leukocytes_ORIGINAL\2011_Nasir';
fileFolder = fullfile(caminho);
dirOutput = dir(fullfile(fileFolder, '*.bmp'));
fileNames = {dirOutput.name}';
numImages = numel(fileNames);

for i = 1 : numImages
    file = [caminho '/' fileNames{i}];
    Im = imread(file);
    

%% LINEAR CONTRAST
R = Im(:,:,1);
G = Im(:,:,2);
B = Im(:,:,3);

R1 = imadjust(R);
B3 = imadjust(B);
G2 = imadjust(G);

F = cat(3,R1,G2,B3);
% imshow(F);
%% HUE COMPONENT
H = rgb2hsv(F);
h = H(:,:,1);

% figure, imshow(h);
%% K MEANS
nrows = (size(h,1));
ncols = (size(h,2));
data = [h(:)];
nColors = 3;
[cluster_idx, cluster_center] = kmeans(data, nColors, 'distance','sqEuclidean',...
  'Replicates',3);
pixel_labels = reshape(cluster_idx, nrows, ncols);
% figure, imshow(pixel_labels,[]), title('image labeled by cluster index');

%% Refazendo a imagem

mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
nuclei_cluster = idx(3);
[u,v] = size(h);
A = 0;
A1 = 0;
for i = 1:u
    for j = 1:v
           if pixel_labels(i,j) == nuclei_cluster
              A(i,j) = 1; 
           else
               A(i,j) = 0;
           end 
    end
end

Med = medfilt2(A, [7 7]);
% [x, y] = getpts(Med);
% reg = regiongrowing(Med,x,y,0.2);
Med = imfill(Med);
% figure, imshow(Med+reg);
for i = 1:u
    for j = 1:v
           if Med(i,j) ~= 1
              R(i,j) = 255;
              G(i,j) = 255;
              B(i,j) = 255;
           end 
    end
end

Final = cat(3,R,G,B);
H1 = rgb2hsv(Final);
%% Segment S band
%%
%%
%%
s = H1(:,:,2);

r = (size(s,1));
c = (size(s,2));
vetor = [s(:)];
C = 3;
[idx, center] = kmeans(vetor, C, 'distance','sqEuclidean','EmptyAction','singleton',...
  'Replicates',3);
label = reshape(idx, r, c);


mean_cluster = mean(center,2);
[tmp, idx] = sort(mean_cluster);
nuclei = idx(3);
[u,v] = size(h);

for i = 1:u
    for j = 1:v
           if label(i,j) == nuclei
              A1(i,j) = 1; 
           else
               A1(i,j) = 0;
           end 
    end
end

% figure, imshow(A1);
    
    
    Med = medfilt2(A1, [7 7]);
    Med = imfill(Med);
    imwrite(Med, file);
end
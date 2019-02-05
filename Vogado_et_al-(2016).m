function [ A ] = WBC_Vogado( f )

transforma = makecform('srgb2lab');
lab = applycform(f,transforma);
L = lab(:,:,1);
A = lab(:,:,2);
B = lab(:,:,3);
AD = imadd(L,B);
f = im2double(f);
r = f(:,:,1);

r = imadjust(r);

g = f(:,:,2);
g = imadjust(g);
b = f(:,:,3);
b = imadjust(b);
c = 1-r;
c = imadjust(c);
m = 1-g;
m = medfilt2(m, [9 9]); % updated in 13/04/2016 - 6:15

m = imadjust(m);
m = imadjust(m);

y = 1-b;
y = imadjust(y);
AD = mat2gray(B);
AD = medfilt2(AD, [9 9]); 
sub = imsubtract(m,AD);
CMY = cat(3,c,m,y);
F = cat(3,r,g,b);
% imshow(sub);
%%
LEN = 21;
THETA = 11;
PSF = fspecial('motion', LEN, THETA);
wnr1 = deconvwnr(sub,PSF,0);
% figure, imshow(wnr1);
%%
ab = double(CMY(:,:,1:3));
nrows = (size(ab,1));
ncols = (size(ab,2));
ab = reshape(ab,nrows*ncols,3);

nrows = (size(c,1));
ncols = (size(c,2));
[x, y] = size(c);
data = [sub(:)];
nColors = 3; % Quantidade de grupos
[cluster_idx, cluster_center] = kmeans(data, nColors, 'distance','sqEuclidean',...
  'Replicates',4);%Aplico o Kmeans
pixel_labels = reshape(cluster_idx, nrows, ncols);
mean_cluster_value = mean(cluster_center,2);
[tmp, idx] = sort(mean_cluster_value);
nuclei_cluster = idx(3);
[u,v] = size(c);
for i = 1:u
    for j = 1:v
           if pixel_labels(i,j) == nuclei_cluster
              A(i,j) = 1; 
           else
              A(i,j) = 0;
           end 
    end
end

se1 = strel('ball',1,1);
se = strel('disk',2);
A = imdilate(A,se1);
A = imerode(A,se);
A = bwareaopen(A, 300);
% figure, imshow(A);
end


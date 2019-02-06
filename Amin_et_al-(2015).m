function [ PP ] = Amin_WBC( Im )
% Im = imread('Im018_1.tif');
R = Im(:,:,1);
%% Pre-processing
HSV = rgb2hsv(Im);
H = HSV(:,:,1);
S = HSV(:,:,2);
V = HSV(:,:,3);
V_equalize = histeq(V);
HSV1 = cat(3,H,S,V_equalize);
%imshow(Im);
RGB = hsv2rgb(HSV1);
%imshow(RGB);
%% Nucleus Segmentation
HSV_NS = rgb2hsv(RGB);
ab = double(HSV_NS(:,:,1:2));
nrows = (size(ab,1));
ncols = (size(ab,2));
ab = reshape(ab, nrows*ncols,2);
nGrupos = 4;
   %% K-means
[idx, center] = kmeans(ab, nGrupos, 'distance', 'sqEuclidean',...
    'Replicates',4);
pixel_labels = reshape(idx, nrows, ncols);
%% Calculo vermelho
i = 1;
soma2 = [1 2 3 4];
for j = 1:nGrupos
      soma(2,j) =  sum(R(pixel_labels == j));
      suma(j) = sum(pixel_labels(pixel_labels == j)); 
end
soma(1,:) = soma2;
soma(3,:) = suma;
[x, y] = size(soma);
for i = 1:x
    for j = 1:y
%         if (soma(i,j) == minimo)
            valor = soma(2,j);
            media(2,j) = valor/soma(3,j);
%         end
    end
end
media(1,:) = soma2;
minimo = min(media(2,:));
[a,b] = size(media);
for i = 1:a
    for j = 1:b
        if (media(i,j) == minimo)
            val = media(1,j);
        end
    end
end

%% Pegando o nucleo
[t,r] = size(pixel_labels);
for i = 1:t
    for j = 1:r
       if(pixel_labels(i,j) == val) 
        Final(i,j) = 0;
       else
        Final(i,j) = 1;
       end
    end
end


mean_cluster_value = mean(center,2);
[tmp, idx] = sort(mean_cluster_value);
nuclei_cluster = idx(3);

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
%figure, imshow(A);
%% Pos-processamento
A = imcomplement(A);
se = strel('disk',7);
PP = imopen(A, se);
PP = imclose(PP, se);
%figure, imshow(PP);

end


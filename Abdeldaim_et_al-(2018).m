function [ BW2 ] = WBC_Abdeldaim( Imagem )
%% Conversão de RGB para CMYK
c = rgb2cmyk(Imagem);

Ciano = c(:,:,1);
Magenta = c(:,:,2);
Amarelo = c(:,:,3);
Preto = c(:,:,4);
%% Histogram equalization thresholding
equalizado = histeq(Amarelo);
ajustado = imadjust(Amarelo);
% figure, imshow(equalizado);

%% Algoritmo de Zack para segmentar a imagem
%I = imcomplement(ajustado);
%figure, imshow(I);

[lehisto x]=imhist(equalizado);
[level]=triangle_th(lehisto,256);
I_bw=im2bw(equalizado,level);
I_bw = imcomplement(I_bw);
% figure, imshow(I_bw);
%% Ajustando a imagem
% pi = 3.14159;
% stats = regionprops('table',I_bw,'Area','ConvexArea','Solidity','Perimeter');
% [u,v] = size(stats);
% v = 1;
% for i = 1:u
%     for j = 1:v
%     Roundness(i,j) = (4*stats.Area(i,j)*pi)/(stats.Perimeter(i,j)^2);
%     end
% end
%Roundness = (4*pi*stats.Area)/(stats.Perimeter);

% Label the blobs.
labeledImage = bwlabel(I_bw);
measurements = regionprops(labeledImage,'Area','Perimeter');
% Do size filtering and roundness filtering.
% Get areas and perimeters of all the regions into single arrays.
allAreas = [measurements.Area];
allPerimeters = [measurements.Perimeter];
% Compute circularities.
Roundness = allPerimeters.^2 ./ (4*pi*allAreas);
Roundness = Roundness';

%% Removendo áreas 

CC = bwconncomp(I_bw);
S = regionprops(CC, 'Area', 'Solidity');
L = labelmatrix(CC);
Roundness = Roundness/10;
Sem_circularidade = ismember(L, find(Roundness > 0.3));
%figure, imshow(Sem_circularidade);

CC = bwconncomp(Sem_circularidade);
S = regionprops(CC, 'Area', 'Solidity');
L = labelmatrix(CC);
Sem_Solidez = ismember(L, find([S.Solidity] >= 0.3));
%figure, imshow(Sem_Solidez);

BW2 = imfill(Sem_Solidez,'holes');
%figure,imshow(I_bw);
%figure,imshow(BW2);


end


function [BW2] = Sedat_et_al(I)
%% Step 0: The original image is taken.
%% Step 1: The image’s intensity values are mapped to a new range by using imadjust.
Im = I;
R = Im(:,:,1);
G = Im(:,:,2);
B = Im(:,:,3);
R = imadjust(R);
B = imadjust(B);
G = imadjust(G);
I = cat(3,R,G,B);
%% Step 2: The RGB image is converted to the grayscale image by using rgb2gray.
RGB = rgb2gray(I);
%% Step 3: The complement of the image is computed by using imcomplement.
gray = imcomplement(RGB);
%% Step 4: Otsu’s method is used to automatically convert
%  the grayscale image to the binary image. The global
%  threshold (level) is computed by using graythresh.
level = graythresh(gray);
BW = im2bw(gray,level);
%% Step 5: The imdilate function dilates the binary image
%  by using the flat, disk-shaped structuring element with
%  radius 1. Thus areas of foreground pixels grow in size
%  while holes within those regions become smaller.
se = strel('disk', 1);
SE = imdilate(BW, se);
% imshow(SE);
%% Step 6:
% gray2 = imcomplement(SE);
gray3 = imcomplement(BW);
BW = imfill(BW, 'holes');
% imshow(BW);

%% Step 7
CC = bwconncomp(BW);
L = labelmatrix(CC);
% figure, imshow(BW);
%% Step 8
stats = regionprops(L, 'Area', 'BoundingBox', 'MajorAxisLength', 'MinorAxisLength', 'PixelList', 'PixelIdxList');

%% Step 9
[u, v] = size(stats);
if (u > 1) 

    Average_major_axis_length = 0;
    Average_minor_axis_length = 0;
    Average_comp = zeros(1,u);
    for i = 1:u
            Average_comp(i) = (stats(i).MajorAxisLength + stats(i).MinorAxisLength)/2;
            Average_major_axis_length = Average_major_axis_length + stats(i).MajorAxisLength;
            Average_minor_axis_length = Average_minor_axis_length + stats(i).MinorAxisLength;
    end

    Average_major_axis_length = Average_major_axis_length / u;
    Average_minor_axis_length = Average_minor_axis_length / u;

    Average_axis_length_for_image = (Average_major_axis_length + Average_minor_axis_length) / 2;

    for i = 1:u
            size_comp(i) = 100*(Average_comp(i)/Average_axis_length_for_image);
    end
    x = (70 * Average_axis_length_for_image)/100;
    %% Step 10
 
cont = 0;
for i = 1:u
    if (size_comp(i) > x)
       BW = bwareaopen(BW, stats(i).Area);
       cont = cont + 1;
    end
end
end

%% Step 11 
se = strel('disk', 1);
BW2 = imerode(BW, se);
% imshow(BW2);
% Z = immultiply(BW2,gray);
end


function [ mask ] = Mohammed_f( FileNames )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% clc
% close all
% clear all
% clear all
% Files=dir('C:\Users\Emad\Dropbox\CLL_BloodImages\oldCLLandNormal');  
                %C:\Users\Mostafa\Dropbox\CLL_BloodImages\CLL images 
                %C:\Users\Mostafa\Dropbox\CLL_BloodImages\oldCLLandNormal
                %C:\Users\Mostafa\Dropbox\CLL_BloodImages\normal lymphocyte images
                %C:\Users\emad\Dropbox\CLL_BloodImages\oldCLLandNormal 
                %C:\Users\Emad\Dropbox\CLL_BloodImages\CLL images number 2\CLL
                %C:\Users\Emad\Dropbox\CLL_BloodImages\CLL images number 2\neutrophils
                %C:\Users\Emad\Dropbox\CLL_BloodImages\CLL images number 2\nomal reactive lymphocytes
                %C:\Users\Emad\Dropbox\CLL_BloodImages\CLL images number 2\normal resting lymphocytes
                %C:\Users\Emad\Dropbox\CLL_BloodImages\CLL images number 2\smudge cells
% [x,y] = meshgrid(-7:7,-7:7);
% sigma=3.5; 
% gxy=exp(-0.5*(x.^2+y.^2)/sigma^2);
% %gxy=ones(7,7);
% gxyn=gxy/sum(gxy(:));

%figure;surf(gxyn)
   
% for k=3:length(Files)
% FileNames=Files(k).name

% FileNames='Im001_1.jpg';

 ColorInputImage=imread(FileNames);
%figure,  imshow(ColorCoinsImage)

GrayImage=rgb2gray(ColorInputImage);%convert the image to gray scale
% figure,  imshow(CoinsImage)
tic
% figure,  imshow(w,[])
BW = im2bw(GrayImage,graythresh(GrayImage));%threshold the image using Otsu's method

BW = edge(BW,'canny');% apply canny edge detector
%  figure, imshow(BW), title('binary gradient mask');%display the binary gradient mask
BW = imdilate(BW,strel('disk',1, 0) );% dilate the binary gradient mask in order to close the open circles
% figure, imshow(BW), title('dilated gradient mask');%display the dilated gradient mask
BW = imfill(BW, 'holes');% fill the circular coins shape
% figure, imshow(BW);title('binary image with filled holes');%binary image with filled holesk
BW = imerode(BW,strel('disk',1,0));% Erode the image to remove the isolated pixels further usage of bwareaopen function maybe used to remove small areas
% figure, imshow(BW), title('final segmented image');%display final segmented image
% BW = bwmorph(BW,'clean');
% w=conv2(double(BW),gxyn,'same');
% BW=uint8(w);

D       = -bwdist(~BW,'chessboard');% 'cityblock', 'chessboard', 'quasi-euclidean', or 'euclidean'
D(~BW)  =-Inf;% min(D(:)); %best value to reduce the effect of local minima
D       = imhmin(D,2,8);% 2% is the height threshold for suppressing shallow 
% D = imimposemin(D,BW); 
L = watershed(D,8);
toc
%L = imfilter(L, ones(1, 1)/1, inf);
%   figure,imshow(label2rgb(L,'jet','k'))
C = size(GrayImage)/2;
s = regionprops(L, 'centroid');
centroids = cat(1, s.Centroid)
% centroids=centroids(3:end,:);
% area = cat(1, s.Area);
% area=area(3:end);
center=size(centroids);
R=C(ones(center(1),1),:);
EquDistance = sqrt((centroids-R).^2) % calculate the Euclidean Distance for every feature in the feature space for the  

%[minValue minIndex] = min(EquDistance)%get the minimum euclidean distance for every features

% [maxarea maxAreaIndex]=max(area)
% if(area(minIndex(1))> area(minIndex(2)))
%      index=minIndex(1);
% else index=minIndex(2);
% end

[dist mindistIndex]=min(EquDistance(2:end,1)+EquDistance(2:end,2))

if isempty(dist)
    mindistIndex=-1;
end
% index;
BWoutline = bwperim(L==mindistIndex+1);%max(minIndex)maxAreaIndex+1
% h=bwareaopen(L,min(area),8);
% figure,imshow(label2rgb(L,'jet','k'))
label=label2rgb(L,'jet','k');
mask              = zeros(size(GrayImage));
mask(BWoutline)   = 1;
mask=imfill(mask,'holes');
%figure;imshow(mat2gray(mask),[]),title('Filled Image');
Segout  = GrayImage;
Segout1 = ColorInputImage;
Segout1(:,:,1)=Segout1(:,:,1).*uint8(mask);
Segout1(:,:,2)=Segout1(:,:,2).*uint8(mask);
Segout1(:,:,3)=Segout1(:,:,3).*uint8(mask);
%figure, imshow(Segout.*uint8(mask)), title('outlined original image');
% %%%%%%%%%%%%%%%%%%
% IC=ColorInputImage;
% I1=rgb2gray(IC);
% I = im2bw(I1,1-graythresh(I1));
% BWs = edge(I,'canny');
% %figure, imshow(BWs), title('binary gradient mask');
% se0  = strel('disk',1, 0);
% BWsdil =imdilate(BWs,se0 );
% % figure, imshow(BWs), title('dilated gradient mask');
% BWdfill = imfill(BWsdil, 'holes');
% %  figure, imshow(BWdfill), title('binary image with filled holes');
% BWnobord = imclearborder(BWdfill,4);
% % figure, imshow(BWnobord), title('cleared border image');
% seD = strel('disk',1,0);
% BWfinal  = imerode(BWnobord,seD);
% seD = strel('disk',1,0);
% BWfinal=imclose(BWfinal,seD);
% imshow(BWfinal);
% %  figure, imshow(BWfinal), title('segmented image');
% ImageProperities  = regionprops(BWfinal, 'ALL');% calculate the area
% Areas     = cat(1, ImageProperities.Area);
% Areamin   = min(Areas);
% Areamax   = max(Areas);
% BWfinal   = bwareaopen(BWfinal,Areamax-Areamin,8);
% BWoutline = bwperim(BWfinal);
% mask      = zeros(size(I1));
% mask(BWoutline)   =1;
% mask=imfill(mask);
% Segout2 = IC;
% Segout2(:,:,1)=Segout2(:,:,1).*uint8(mask);
% Segout2(:,:,2)=Segout2(:,:,2).*uint8(mask);
% Segout2(:,:,3)=Segout2(:,:,3).*uint8(mask);
%%%%%%%%%%%%%%%%%%
% im=[ColorInputImage Segout1;Segout2 Segout1-Segout2];
%  figure, imshow(im,[])% instead of label you can use (Segout1-Segout2)
%  title(sprintf(' %10s outlined original image ',FileNames))%title('outlined original image');
% j='P.jpg'; 
% figure, imshow(Segout1);
% name=[FileNames(:,1:end-3),j];
% imwrite(uint8(im),name);
%k-2
% kk='GT';
% nname=[kk,FileNames];
% im2=imread(nname);
% Inuc=(IC-im2);
% It=Inuc;
% %figure,imshow(Inuc,[]),title('InucSeg')
% %figure,imshow(Segout1,[]),title('Segout1')
% error=(Inuc-Segout2);
% %figure,imshow(error,[]),title('error')
% Inuc=rgb2gray(Inuc);
% Inuc = im2bw(Inuc,graythresh(Inuc));
% %figure,imshow(double(Inuc),[]),title('bw')
% InucTot=sum(Inuc(:));
% error=rgb2gray(error);
% error = im2bw(error,graythresh(error));
% max(error(:));
% max(Inuc(:));
% %figure,imshow(double(error),[]),title('ebw')
% errorTot=sum(error(:));
% AccuracyNuc(k-2,:)=(1-(errorTot/InucTot))*100;
% II=im2-It;
% % figure,imshow(II,[]),title('Nuc')
% %%%%%%%%%%%%%%%%%%%%%
% kk='GTC';
% nname=[kk,FileNames];
% im3=imread(nname);
% Inuc=(IC-im3);
% Icell=Inuc;
% %figure,imshow(Inuc,[]),title('InucSeg')
% %figure,imshow(Segout1,[]),title('Segout1')
% error=(Inuc-Segout1);
% %figure,imshow(error,[]),title('error')
% Inuc=rgb2gray(Inuc);
% Inuc = im2bw(Inuc,graythresh(Inuc));
% %figure,imshow(double(Inuc),[]),title('bw')
% InucTot=sum(Inuc(:));
% error=rgb2gray(error);
% error = im2bw(error,graythresh(error));
% max(error(:));
% max(Inuc(:));
% % figure,imshow(double(error),[]),title('ebw')
% errorTot=sum(error(:));
%  AccuracyCell(k-2,:)=(1-(errorTot/InucTot))*100;
% % all(k-2,:)=Files(k).name;
% II=im2-Icell;
% %figure,imshow(II,[]),title('Cyto')
% %%%%%%%%%%%%%%%%%%%%%
% %  Icyto=Icell-It;
% Icytoseg=Segout1-Segout2;
% %  Idiff=Icyto-Icytoseg;
% 
% Icyto=Segout1-It;
% Idiff=Icyto-Icytoseg;
% 
% Idiff=rgb2gray(Idiff);
% Icyto=rgb2gray(Icyto);
% Idiff = im2bw(Idiff,graythresh(Idiff));
% Icyto = im2bw(Icyto,graythresh(Icyto));
% Idifft=sum(Idiff(:));
% Icytot=sum(Icyto(:));
% AccuracyCyto(k-2,:)=(1-(Idifft/Icytot))*100;
% %figure,imshow(Idiff,[])
%  end

end


function out_mask = WBC_SegMadhloom(input_image, show_images)
%This function is built to segment the nuclei of white blood cells
% the input is an image and an option to show the resluts on screen or not (1 or 0)
%and the number of white blood cells
% For more information abou the algorithm please refere to: 
% "An Efficient Technique for White Blood Cells Nuclei Automatic Segmentation", 
% IEEE SMC conference 2012
%  The algorithm in this function is an implementation if the paper:
% H. T. Madhloom, S. A. Kareem, H. Ariffin, A. A. Zaidan, H. O. Alanazi and 
% B. B. Zaidan, "An Automated White Blood Cell Nucleus Localization and 
% Segmentation using Image Arithmetic and Automatic Threshold," Journal of 
% Applied Sciences, vol. 10, no. 11, pp. 959-966, 2010. 
if(ischar(input_image)==0)
    a1 = input_image;
     input_image = 'Memory';
else
    a1=imread(input_image);
end;
if show_images == 1, figure;imshow(a1);title('Input image'); end
%% Convert to gray scale image
a=rgb2gray(a1);
if show_images == 1, figure;imshow(a);title('Input image in gray scale'); end;
%% Adjust image intensity values with a linear contrast stretching 
L=imadjust(a);
if show_images == 1, figure;imshow(L);title('L: Input image with intensity values adjusted by a linear contrast stretching'); end;
%% Enhance contrast using histogram equalization 
H = histeq(a);
if show_images == 1, figure;imshow(H);title('H: Enhance contrast using histogram equalization'); end
%% Brighten most of the details in the image except the nucleus
R1=imadd(L,H);
if show_images == 1, figure;imshow(R1);title('R1: Brighten most of the details in the image except the nucleus'); end
%% Highlight all the objects and its borders in the image including the cell nucleus
R2 = imsubtract(R1,H);
if show_images == 1, figure;imshow(R2);title('R2: Highlight all the objects and its borders in the image including the cell nucleus'); end
%% Remove almost all the other blood components while retaining the nucleus with minimum affect of distortion on the nucleus part of the white blood cell.
R3=imadd(R1,R2);
if show_images == 1, figure;imshow(R3);title('R3: Remove almost all the other blood components while retaining the nucleus with minimum affect of distortion on the nucleus part of the white blood cell.'); end
%check histogram
% figure;
% hi=imhist(R3);
% hi1=hi(1:2:256);
% horz=1:2:256;
% bar(horz,hi1);
% % %axis([0 255 0 1400]);
% % %set(gca, 'xtick', 0:50:255);
% % %set(gca, 'ytick', 0:2000:15000);
%  xlabel('Gray level' );
% ylabel('No of pixels' );
% title('Histogram before opening the image');

% 'Time elapsed to calculate image preparations (algorithm)' 
% toc(programStart)
filterStart=tic;
%=====================
%implements a 3-by-3 minimum filter
%reduce noise, preserve edges and increase the darkness of the nuclei
for i=1:1
    R3=ordfilt2(R3,1,ones(3,3));
end
% 'Time elapsed to calculate 3-by-3 minimum filter' 
% toc(filterStart)
%figure; imshow(R3);
%check histogram afteraplying the minimum filter
% figure;
% hi=imhist(R3);
% hi1=hi(1:2:256);
% horz=1:2:256;
% bar(horz,hi1);
% % %axis([0 255 0 10000]);
% % %set(gca, 'xtick', 0:50:255);
% % %set(gca, 'ytick', 0:2000:15000);
%  xlabel('Gray level' );
% ylabel('No of pixels' );
% title('Histogram before opening the image');
%=====================================================
%% Global threshold using Otsu’s method 
level=graythresh(R3);
bw=im2bw(R3,level);
%Complement image
bw = imcomplement(bw);
if show_images == 1, figure;imshow(bw);title('Global threshold using Otsu’s method '); end
%% Show the results
out_mask = bw;
if show_images == 1
    %Re-count objects in the image
    cells = bwconncomp(bw);
    no_of_WBCs=cells.NumObjects;
   
    a2 = a1;
    for j=1:no_of_WBCs
        %return the coordinates for the pixels in object j
        [r, c] = find(bwlabel(bw)==j);
        rc = [r c];
        %marking the location of the nuclei of white blood cells by a white
        %colour
        for i=1:max(size(rc))
            a2(rc(i,1),rc(i,2),:)=uint8(a1(rc(i,1),rc(i,2),:)*1.5);
        end
    end
    h = figure(10);
    set(h, 'Name', 'WBC segmentation results','NumberTitle','off', 'OuterPosition',[1 1 1600 600])
    subplot(1,3,1);imshow(a1);
    [~, name ext] = fileparts(input_image);
     title(['The original image [' name ext ']']);

    subplot(1,3,2);imshow(a2);
    title('The original image with the white spot(s) over the WBC nucleus(ei)');
    %show the final segmented image
    bw1 = imcomplement(bw);
    amask = a;
    amask(bw1) = 255;
    subplot(1,3,3);imshow(amask);
    title('The final segmented image mask');
%     disp(['Total program time is: ' num2str(toc(programStart))]);

%     Check Histogram of masked gray scale image
%     figure(5);  % imshow(a(bw))       
%     hi=imhist(amask); % imshow(amask)
%     hi1=hi(1:2:256);
%     horz=1:2:256;
%     bar(horz,hi1);
% %     axis([0 255 0 1400]);
%     % %set(gca, 'xtick', 0:50:255);
%     % %set(gca, 'ytick', 0:2000:15000);
%     xlabel('Gray level' );
%     ylabel('No of pixels' );

%     Save images
    path = fileparts(input_image);
    imwrite(a1, [path '\a1Mah.bmp'], 'bmp')
    imwrite(a2, [path '\a2Mah.bmp'], 'bmp')
    imwrite(bw1, [path '\abwMah.bmp'], 'bmp')
end

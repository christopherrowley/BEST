function [ recombine_seg ] = cerebellum_extract2(current_seg, norm_img2,sag_slices,threshold)
%  where norm_img_masked is the normalized image with an initial cerebrum
%  mask applied to remove all/most skull. 
%  the goal of this script is to erode the center line a bit extra to
%  remove the central dura

%  current_seg = best segmentation up to this point
%  norm_img2 = normalized image 
%  sag_slices = the window size to erode in the center of the image
%  threshold = the GM threshold found from sharp_csf 

%%%%%%%%%%%%%%%%%%%%%%%%
%crop images to take center of image.

% sag_slices = 50;
% norm_img_masked = DS_2.*norm_img;
% current_seg = DS_2;

norm_img_masked = current_seg.* norm_img2;
sag_slices = round(sag_slices/2);

dim_MRimage = size(norm_img_masked);
cropdim =  round(0.5*dim_MRimage);

cropnorm = norm_img_masked(cropdim(1)-sag_slices:cropdim(1)+sag_slices,:,:);

% figure;
% title('upsampled and gaussian blurred - sharpened')
% imshow3Dfull(cropnorm)

%find new threshold to remove midline:
threshold3 = threshold;

crop_mask = zeros(size(cropnorm));
crop_mask(cropnorm > threshold3) = 1;

% figure;
% title('further eroded midline check')
% imshow3Dfull( (crop_mask+1) .*cropnorm)

CC = bwconncomp(crop_mask,6);
L = labelmatrix(CC);
%want just the brain, calculate areas, we want the one with largest area
S = regionprops(CC, 'Area');
idx = find([S.Area] == max([S.Area])); %find the largest area value
L3 = zeros(size(L));
L3(L == idx) = 1; 

% se = strel3d(3); 
% dilate_seg = imdilate(L3,se);

recombine_seg = current_seg; 
recombine_seg(cropdim(1)-sag_slices:cropdim(1)+sag_slices,:,:) = L3; %dilate_seg;

% figure;
% title('further eroded midline check')
% imshow3Dfull( (L3+1) .*cropnorm)


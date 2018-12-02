function [g_blur2] = rough_skull_strip(norm_img)
%% Try to generate an inital skull segmentation from the image. 

% %using norm_image
% figure;
% imshow3Dfull(norm_img)

[x,y,z] = size(norm_img);
midx = round(x/2); %moving from left to right
midy = round(y/2); %move from back to front of head
midz = round(z/2); %move from neck to top of head

%binarize based on WM intensity 
%take a cropped central region to assess due to background noise
wm_crop = zeros(x,y,z);
wm_crop(midx-40:midx+40,midy-40:midy+40,midz-40:midz+40) = norm_img(midx-40:midx+40,midy-40:midy+40,midz-40:midz+40); 
wmintensity = threshold_finder_gaussfit(wm_crop,'wm_peak') - 0.1;

%wmintensity = 0.88;
norm_img2 = norm_img;
norm_img2(norm_img > 1.4) = 0; %this should help prevent skull from being largest object. was 1.4

bin_mri = zeros(x,y,z);

if wmintensity > 1
    wmintensity = 1;
end

for i = 1:z
    bin_mri(:,:,i) = im2bw(norm_img2(:,:,i), wmintensity); %switch to imbinarize on newer matlab
end

bin_mri2= imfill(bin_mri,'holes');

% figure;
% imshow3Dfull(bin_mri2)

%create a box to remove the brain stem area so it isn't connected to
%anything. 
rm_stem = ones(x,y,z);
rm_stem(midx-15:midx+15,midy-30:midy+20,1:midz-10) = 0; %brain stem
rm_stem(midx-20:midx+20,midy-45:midy+5,1:midz-25) = 0; %part of cerebellum
rm_stem(midx-20:midx+20,midy-80:end,1:midz-45) = 0; % inferior to orbital frontal cortex
rm_stem(1:end,1:midy-10,1:midz-45) = 0; % inferior to orbital frontal cortex


% 
% figure;
% imshow3Dfullseg(norm_img,[0 2],rm_stem)

bin_mri_no_stem = bin_mri2 .* rm_stem;
bin_mri_no_stem = bwareaopen(bin_mri_no_stem, 1000,6);
% 
% figure;
% imshow3Dfull(bin_mri_no_stem)


% brain box
bb = zeros(x,y,z);
bb(midx-15:midx+15,midy-30:midy+30,midz:midz+25) = 1; %brain stem

CC = bwconncomp((bin_mri_no_stem),6); % calculates the connected components
L = labelmatrix(CC); % creates an image where all CC are one value
idx = max(max(max(double(L) .* double(bb)))); %find value of WM component
L2 = zeros(size(L));
L2(L == idx) = 1; 

% 
% figure;
% imshow3Dfullseg(bigObjects_erode,[0 1],L2)

%do a large dilate to be sure to add in all brain! 
se = strel3d(4); %might work better as a 2 step, in order to better get gyral crowns %was 8
d1 = imdilate(L2,se);

se = strel3d(5); % was 6 
d2 = imdilate(d1,se);

se = strel3d(5); % was 6 
d3 = imdilate(d2,se);

 g_blur = double(imgaussfilt3(d3,5)); 
 g_blur2 = round (g_blur);


%Gaussian blur and then round should perform a dilation,while adding in
%anything removed that was between cortex. 

% figure;
% imshow3Dfull((1+g_blur2).*norm_img,[0 4])






















function [ mask_minus_sinus ] = sag_sinus_rm(current_seg, norm_img2)

%try a convolution method with 3 filters? 
%troubleshooting
% norm_img2 = norm_img;
% current_seg = L2;

% load the 3 kernels for image convolution, need odd dimensions 
filename = '/home/chris/Desktop/segmentation_prg/conv_kern.xlsx';

tri_kern = xlsread(filename,'triangle_closed');
circ_kern = xlsread(filename,'circle_closed');

tri_kern = imgaussfilt(tri_kern);
circ_kern = imgaussfilt(circ_kern);

% figure;
% imagesc(tri_kern)
% axis equal
% figure;
% imagesc(circ_kern)
% axis image

% perform cross correlation on axial slice
% crop the image to midline to speed up and reduce chance of false positives
cropdim =  round(0.5*(size(current_seg)));
sag_slices = 18;
cropnorm_sag = current_seg(cropdim(1)-sag_slices:cropdim(1)+sag_slices,:,:);
cropnorm_sag = imgaussfilt3(cropnorm_sag);

% figure;
% imshow3Dfull(cropnorm_sag)

%only run if enough white pixels (to prevent wasted time on empty slices)
kern_tri_size = size(tri_kern);
kern_circ_size = size(circ_kern);
cropsag_dim = size(cropnorm_sag);

cross_tri = zeros(cropsag_dim(1) + (kern_tri_size(1) -1), cropsag_dim(2) + (kern_tri_size(2) -1), cropsag_dim(3));
cross_circ = zeros(cropsag_dim(1) + (kern_circ_size(1) -1), cropsag_dim(2) + (kern_circ_size(2) -1), cropsag_dim(3));

for i = 1 :size(cropnorm_sag, 3) 
    if sum(sum(cropnorm_sag(:,:,i))) >10
            cross_tri(:,:,i) = normxcorr2new(tri_kern,cropnorm_sag(:,:,i));
            cross_circ(:,:,i) =normxcorr2new(circ_kern,cropnorm_sag(:,:,i));
    end % if not enough non-zero pixels, then it remains at 0
end

% figure;
% imshow3Dfull(cross_tri)
% figure;
% imshow3Dfull(cross_circ)
% figure;
% imshow3Dfull(cross_rect)

cross_tri_thres = zeros(size(cross_tri));
cross_tri_thres(cross_tri >= 0.5) = 1;
cross_circ_thres = zeros(size(cross_circ));
cross_circ_thres(cross_circ >= 0.5) = 1;

% figure;
% imshow3Dfull(cross_tri_thres)
% figure;
% imshow3Dfull(cross_circ_thres)
% % figure;
% % imshow3Dfull(cross_rect_thres)

resized_tri = post_xc_resize(tri_kern, cross_tri_thres);
resized_circ = post_xc_resize(circ_kern, cross_circ_thres);
% resized_rect = post_xc_resize(rect_kern, cross_rect_thres);

combo_info = resized_tri + resized_circ;

% figure;
% imshow3Dfull(combo_info)

combo_info(combo_info < 2) = 0;
combo_info(combo_info>1) = 1;

% figure;
% imshow3Dfull(combo_info)

se = strel('rectangle',[15 15]); 
dil_comb_info = zeros(size(combo_info));
for i = 1 :size(combo_info, 1) 
    if sum(sum(combo_info(i,:,:))) >1
        dil_comb_info(i,:,:) = imdilate(combo_info(i,:,:),se);
    end % if not enough non-zero pixels, then it remains at 0
end

dil_comb_info2 = imgaussfilt3(dil_comb_info,3);
dil_comb_info2(dil_comb_info2 >0.2) = 1; 
dil_comb_info2(dil_comb_info2 <= 0.2) = 0; 
% figure;
% imshow3Dfull(dil_comb_info2)

CC = bwconncomp(dil_comb_info2,6);
%want just the brain, calculate areas, we want the one with largest area
S = regionprops(CC, 'Area','PixelIdxList');  %does this always return largest one first? 
[~, idx] = max([S.Area]);
seg_rm = zeros(size(dil_comb_info2));
seg_rm(S(idx).PixelIdxList) = 1; 

% figure;
% imshow3Dfull(seg_rm)

brain_seg = zeros(size(norm_img2));
brain_seg(cropdim(1)-sag_slices:cropdim(1)+sag_slices,:,:) = seg_rm;
se = strel3d(3); 
brain_seg = imdilate(brain_seg,se);


brain_seg2 = brain_seg .* current_seg;
brain_seg2 = bwareaopen(brain_seg2, 20,6); %remove small pieces that may be GM

mask_minus_sinus = current_seg - brain_seg2; 
























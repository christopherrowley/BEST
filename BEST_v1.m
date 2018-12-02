%Main Script
% currently requires to be registered to MNI152 space, then fslorient
% -deleteorient run on it. 
clear all
%%%%%%%%%%%%%%%%%%%%%%
tic

SUBJECT = {'CH24_20141104'}; %      }; 
% ,'May_8_2015','May_29_2015','Me_December_20_2016'}; %'image1' 'image2' 'image4'
% for exercise study:  'JO23_20140703' 'JU19_20140701' 'MI30_20140701'
% 'May_29_2015' missing orbital frontal cortex and some at the top of the
% brain, potentially being cropped out in early stages. 

DATADIR = '/media/chris/6096EFEF96EFC39E/user/Research/Aerobic_Exercise_MRI/data/';
%T1filename =  strcat(DATADIR,'image3/T1image.nii');
%DATADIR = '/media/chris/6096EFEF96EFC39E/user/250micron_data/';
%DATADIR = '/media/chris/6096EFEF96EFC39E/user/Research/Aerobic_Exercise_MRI/trial/';

%pause_time = 100; %time you need to check segmentations before it runs next sub. 
%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%
for z = 1:size(SUBJECT,2)
    
T1filename =  strcat(DATADIR,SUBJECT{z},'/image/Synth_T1_orient.nii'); %register/Synth_T1_orient.nii' %'/bravo_standard.nii');
left_filename = strcat(DATADIR,SUBJECT{z},'/image/BEST_left_gm.nii');
right_filename = strcat(DATADIR,SUBJECT{z},'/image/BEST_right_gm.nii');
left_wmfilename = strcat(DATADIR,SUBJECT{z},'/image/BEST_left_wm.nii');
right_wmfilename = strcat(DATADIR,SUBJECT{z},'/image/BEST_right_wm.nii');

%final_filename = strcat(DATADIR,SUBJECT{z},'/register/norm_img.nii');

%%%%%%%%%%%CODE%%%%%%%%%%%

%load image
T1image = load_nii(T1filename); 
T1image.img = double(T1image.img);

%FIX IMAGE LIMITS
T1image.img(isnan(T1image.img)) = 0 ; 
T1image.img(isinf(T1image.img)) = 0 ;%remove nan's 

%reorient image if not down previously
%T1image.img = flip(T1image.img,2); %BD data with delete orient needs this
%T1image.img = flip(T1image.img,3);
% figure;
% imshow3Dfull(T1image.img )

%%%%%%%%%%%%%%%%%%%%%%% Get the intensity self-corrected/normalized image
norm_img = PV_divison(T1image.img); %intensity correction before downsample

figure;
imshow3Dfull(norm_img,[0,2])
% 
% %%%%%%%%%%%Next remove unneeded information. NOT NECESSARY WITH PV_DIVSION FIX 
% norm_img(1:18,:,:) = 0;
% norm_img(165:182,:,:) = 0;
% norm_img(:,1:20,:) = 0;
% norm_img(:,200:218,:) = 0;

init_seg = rough_skull_strip(norm_img);
mask_minus_sinus = sag_sinus_rm_binning(init_seg, norm_img);
MRimage = norm_img .* mask_minus_sinus;

% figure;
% imshow3Dfullseg(norm_img,[0 2],init_seg)
% figure;
% imshow3Dfullseg(norm_img,[0 2],mask_minus_sinus)
% figure;
% imshow3Dfull(MRimage)

%upsample image by 2
dim_MRimage = size(MRimage);
UPS_dim = dim_MRimage.*2; 
UPSAMP_2 = resize(MRimage,UPS_dim); 

% figure;
% title('upsampled')
% imshow3Dfull(UPSAMP_2)

%%image sharpening
%everything has to sum to 0
filt2(:,:,1) = [-0.5 -0.5 -0.5; -0.5 -0.5 -0.5; -0.5 -0.5 -0.5]; % 4.5
filt2(:,:,2) = [-0.5 -0.5 -0.5; -0.5 13 -0.5; -0.5 -0.5 -0.5]; %1.5;1 x; 1.5
filt2(:,:,3) = [-0.5 -0.5 -0.5; -0.5 -0.5 -0.5; -0.5 -0.5 -0.5]; %4.5 

sig2 = double(imgaussfilt3(UPSAMP_2,1));
highpassfilt = convn(sig2,filt2);
highpassfilt = highpassfilt(1:end-2,1:end-2,1:end-2); %convolution adds two rows, so remove them

% figure;
% title('highpass filtered')
% imshow3Dfull(highpassfilt)

img_sharp = highpassfilt + UPSAMP_2; %sharpened image is the normal image plus a highpass filtered image

%%% %%% %%% %%%  remove files to reduce RAM impact
% clear highpassfilt sig2 MRimage mask_minus_sinus init_seg
%%% %%% %%% %%% 

%sharpened image gives a better separation of cortex from dura
 %%%% Sharp image gives a steep inflection point between CSF and GM. find
 %%%% that and threshold on it!!! 
 
threshold = threshold_finder(img_sharp,'sharp_csf');  %want to erode a bit inwards, so increase threshold.  %used to be +0.35, but now threshold_finder shifts it automatically towards wm peak
sharp_seg = zeros(size(img_sharp));
sharp_seg(img_sharp > threshold) = 1; %remove based on lower threshold

% figure;
% imshow3Dfullseg(UPSAMP_2, [0 2], sharp_seg)

%we want a thresholded norm image to make final decisions on what to keep
threshold1 = threshold_finder(norm_img,'norm_csf');  %was UPSAMP_2 but then if could be too low. 

%%%%%%%%%%%%%%%%%%%%%%
%% Temporary Fix GM threshold
%%%%%%%%%%%%%%%%%%%%%%
% threshold1 = 0.48;

norm_seg = zeros(size(UPSAMP_2));
norm_seg(UPSAMP_2 > threshold1) = 1; %remove based on lower threshold

% figure;
% imshow3Dfullseg(UPSAMP_2, [0 2], norm_seg)

%First erode sharp_seg to remove floating dura
%%%% Need to add in WM as holes can come up from sharpened image which
%%%% creates issues...
sharp_seg2 = imfill(sharp_seg,'holes');

% figure;
% imshow3Dfull(sharp_seg2)

se = strel3d(5); %use circular to keep round
eroded_img_seg = imerode(sharp_seg2, se); %erode image
bigObjects = bwareaopen(eroded_img_seg, 10000,6); %remove small objects likely to be dura
% 
% figure;
% imshow3Dfull(bigObjects)

%second dilate
se = strel3d(5); 
dilate_skull_seg = imdilate(bigObjects,se); %apply dilation to replace cortex. 

newpix = dilate_skull_seg .* norm_seg; %only keep new pixels from dilation if we had them before.  

%%% %%% %%% %%% remove files to reduce RAM impact
% clear sharp_seg sharp_seg2 norm_seg bigObjects dilate_skull_seg eroded_img_seg img_sharp
%%% %%% %%% %%% 
smoothed_seg = ordfilt3(newpix,'med',5); %median filter. Looks beautiful

%may be underestimated , so fill in cortex based on edges, then threshold. 
filled_c =fill_cortex(smoothed_seg);  

temp = filled_c .* UPSAMP_2;
temp(UPSAMP_2 < threshold1) = 0; %remove CSF added from smoothing
final_ups_mask = zeros(size(temp));
final_ups_mask(temp>0) =1;

%%
%upsample contour mask. 
DS_2 = resize(final_ups_mask, size(norm_img)); %return the image to its proper size.
DS_2 = round(DS_2);
% DS_3 = DS_2;

DS_3 = cerebellum_extract(DS_2, norm_img); %still having issues. Likely
%need to code an accurate cerebellum segmentation. 

%DS_4 = cerebellum_extract2(DS_3, norm_img ,50, threshold1+0.1);

% figure;
% imshow3Dfullseg(norm_img,[0 2],DS_2);

se = strel3d(3); 
dilate_DS_3 = imdilate(DS_3,se); %some GM missing, so restore. 
dilate_DS_3(norm_img < threshold1) = 0; % remove CSF gain from dilation

smoothed_seg2 = ordfilt3(dilate_DS_3,'med',3); % use contraint that cortex must be smooth
smoothed_seg2(norm_img < threshold1) = 0; %final threshold

%%% %%% %%% %%% remove files to reduce RAM impact
% clear UPSAMP_2 smoothed_seg filled_c final_ups_mask DS_2 DS_3 new_pix temp
%%% %%% %%% %%% 

CC = bwconncomp(smoothed_seg2,6);
L = labelmatrix(CC);
%want just the brain, calculate areas, we want the one with largest area
S = regionprops(CC, 'Area');
idx = find([S.Area] == max([S.Area])); %find the largest area value
L2 = zeros(size(L));
L2(L == idx) = 1; 

figure;
imshow3Dfullseg(norm_img,[0 2],L2);

%%
%now segment the WM using FANTASM. This will also be useful for final GM
%cleaning

wm_seg = FANTASM(L2, norm_img,  100, 0.01, 0.2, 2, 0.05); %inital wm segmentation

[adj_img] = seg_int_corr(L2,wm_seg,norm_img); %use inital segmentations to correct intensity of image

%%% Need to better fix GM intensity near Motor cortex - too bright, WM isnt
%%% good there
% wm_seg2 = FANTASM(L2, adj_img,  100, 0.01, 0.2, 2, 0.15); %recalculate wm segmentation with corrected image
% 
% se = strel3d(3);
% L2erode = imerode(L2, se);
% wm_seg2 = wm_seg2 .* L2erode;

% figure;
% imshow3Dfullseg(adj_img,[0 2],wm_seg2);


%%
% now split brain left and right :)
 %maybe take into account that the central dura is close to iso-intense to
 %cortex. a way to get this? 
 
[left_mask, right_mask] = brain_split4( adj_img, L2);  %%%%% Use Me_December_2016 to trouble shoot a non-straight dividing line between hemispheres

left_mask =double( bwareaopen(left_mask, 10000,6));
right_mask = double( bwareaopen(right_mask, 2500,6));

% figure;
% imshow3Dfullseg(norm_img,[0 2],right_mask);

%%
%with an accurate mask now generated. give it one more shot of growing WM
left_upd = left_mask; 
left_upd(norm_img<threshold1) = 0;
right_upd = right_mask;
right_upd(norm_img<threshold1) = 0;

figure;
imshow3Dfullseg(norm_img,[0 2],right_upd);
figure;
imshow3Dfullseg(norm_img,[0 2],left_upd);

drawnow  %might work better than a pause

%%
%for some reason, it is flipped in both X and Y

% %save Left gm
% savethis = flip(left_upd,2);
% savethis1 = flip(savethis,1);
% savethis2 = make_nii(savethis1);
% save_nii(savethis2,left_filename)
% 
% %save right gm
% savethis = flip(right_upd,2);
% savethis1 = flip(savethis,1);
% savethis2 = make_nii(savethis1);
% save_nii(savethis2,right_filename)
% 
% % %save Left wm
% % savethis = flip(wm_seg2.*left_upd,2);
% % savethis1 = flip(savethis,1);
% % savethis2 = make_nii(savethis1);
% % save_nii(savethis2,left_wmfilename)
% % 
% % %save right wm
% % savethis = flip(wm_seg2.*right_upd,2);
% % savethis1 = flip(savethis,1);
% % savethis2 = make_nii(savethis1);
% % save_nii(savethis2,right_wmfilename)

done = SUBJECT{z} %display who is finished

toc
%pause 
end





























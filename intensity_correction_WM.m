%%%%%%%%%
% this code is used to normalize signal intensity based on the WM
% it runs under the assumption the smoothly varying changes in WM intensity
% are caused by RF shading. This shading will be estimated and extrapolated
% for the nearby GM. 
%%%%%%%%%

% Psuedo Code: 
% 1. erode wm
% 2. Calculate the mean intensity of the eroded wm (used to normalize later)
% 3. calculate distance map
% 4. For i=1:(max(distance_map))
% 	calculate moving window max
% 	new_voxels = zeroes (size(distance_map))
% 	new_voxels(distance_map==i)=1
% 	wm_extend_intensity(new_voxels ==1) = moving_max
% 	end
% 5. Smooth 3d gaussian. 
% 6. wm_extend_intensity /eroded_mean_intensity (should then be centered on 1 relative to the original wm which will provide a proportional factor to adjust the images)
% 
% 


%load image to correct, and WM segmentation. 
SUBJECT = {'005-322'};
DATADIR = '/Users/christopherrowley/Desktop/track-hd/';

for z = 1:size(SUBJECT,2)
    
    T1filename = strcat(DATADIR,SUBJECT{z},'/T1_orient.nii');
    T2filename = strcat(DATADIR,SUBJECT{z},'/T2_orient.nii');
    left_wm_file = strcat(DATADIR,SUBJECT{z},'/left_wm.nii'); %step4-jist/
    right_wm_file = strcat(DATADIR,SUBJECT{z},'/right_wm.nii'); %step4-jist/
    corrected_image_file = strcat(DATADIR,SUBJECT{z},'/T1_T2_intcor.nii');
    field_estimate_file = strcat(DATADIR,SUBJECT{z},'/field_est.nii');

    T1image = load_nii(T1filename); 
    T1image.img = double(T1image.img);
    T2image = load_nii(T2filename); 
    T2image.img = double(T2image.img);
    left_wm_seg = load_nii(left_wm_file); 
    left_wm_seg.img = double(left_wm_seg.img);
    right_wm_seg = load_nii(right_wm_file); 
    right_wm_seg.img = double(right_wm_seg.img);
    
    
    wm_seg = left_wm_seg.img + right_wm_seg.img; %hopefully these are the same size so this works
    T2image.img(T2image.img<0.001) = 0.01;
    t1t2img = T1image.img ./ T2image.img; %assuming gaussian smoothness, dividing the two should leave gaussian shading still. 
    t1t2img(t1t2img<0) = 0 ;
    t1t2img(t1t2img>5) = 5 ;
      
figure;
imshow3Dfull(wm_seg);    
    
figure;
imshow3Dfull(t1t2img);
    
    
    %erode wm
    se = strel(ones(2,2,2)); %check that this isnt too much
    wm_erode = imerode(wm_seg, se); %erode to avoid partial voluming issues
     
figure;
imshow3Dfull(wm_erode);

    %calculate distance map 
    dist_from_wm = ceil(bwdist(wm_erode)); %we are going to move stepwise, so dont need precision
    dist_from_wm = double(dist_from_wm);
    
figure;
imshow3Dfull(dist_from_wm);
    
    %initiate extended wm value 
    wm_extend_intensity = wm_erode .* t1t2img; %start with only the intensity in the eroded WM. This will change iteravely
    final_max_img = wm_extend_intensity;
    
tic
    for i = 1:45  %save time by not going all the way out %max(max(max(dist_from_wm)))
       
        moving_max = moving_window_max(final_max_img, dist_from_wm,i);
        new_voxels = zeros(size(dist_from_wm)); %reset matrix
        new_voxels(dist_from_wm==i)=1; %select the voxels the next step out  
        new_intens = new_voxels .* moving_max; 
        final_max_img = final_max_img + new_intens; %iteratively build up

    end
toc

% 
% figure;
% imshow3Dfull(final_max_img);
% 
% figure;
% imshow3Dfull( new_voxels);    
%   
% figure;
% imshow3Dfull(new_intens);
    
    %becasue we are assuming guassian smoothness in the shading, smooth
    %away
    Iblur = imgaussfilt3(final_max_img,2); %may need to change? it needs to be big enough that it doesnt remove anatomical variation. 

% figure;
% imshow3Dfull(Iblur);
    
    mean_wm = mean(mean(mean(t1t2img(wm_erode==1)))); %normalize based on mean wm intensity
    
    shading_field = Iblur / mean_wm;
    shading_field(shading_field <0.01) =0.1;
    
    corrected_image = t1t2img ./ shading_field; 
    
figure;
imshow3Dfull(shading_field);

figure;
imshow3Dfull(corrected_image);

    %save outputs 
    final_save = T1image; 
    final_save.img = corrected_image;
    save_nii(final_save,corrected_image_file)
    
    final_save.img = shading_field;
    save_nii(final_save,field_estimate_file)

end 






%check the effect of dividing 2 gaussians
% 
% x = -3:.1:3;
% norm1 = normpdf(x,-0.8,0.5); %x, mean, std
% norm2 = normpdf(x,-0.8,0.5);
% norm3 = normpdf(x,1,0.5); %x, mean, std
% norm4 = normpdf(x,0.9,0.5);
% 
% g1 = norm1 + norm3;
% g2 = norm2+norm4;
% 
% div_gaus = g1 ./g2;
% 
% figure;
% plot(x,g1)
% hold on
% plot(x,g2)
% plot(x,div_gaus)











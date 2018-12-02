function [ left_mask, right_mask ] = brain_split4( input_image, finalseg)
%Goal is to define a plane near the center. There should be more 1's in the
%center of the inverse brain mask. Use those as data points to drive a
%2D plane fit. 

% input_image = norm_img;
% finalseg = L2;

T1seg = input_image .* finalseg;

%now in the 3D case, slice by slice

img_size = size(T1seg);
half_size = ceil(img_size/2);
st_idx = half_size(1)-10; %define a narrow search window. 
end_idx = half_size(1)+10;
mp = T1seg(st_idx:end_idx,:,:);

half_brain = zeros(size(finalseg));
stored_vals = zeros(size(finalseg)); %bluring seems to ruin parts with low numbers

for i = 1:size(finalseg,3) %for each slice
%for i = 90:100 %for each slice
    if (sum(sum(mp(1:10,:,i))) >50) && (sum(sum(mp(11:21,:,i))) >50)  %only run the contour if there are lots of datapoints on each side!
        %%%%%%%%%2D
        mp2 = mp(:,:,i);
        low = zeros(size(mp2));

        for x = 1:size(mp2,2)
            num = min(mp2(:,x));
            temp = low(:,x);
            temp(mp2(:,x) ==num) =1; 
            low(:,x) = temp;
        end
        %for each column
        %if one min value exists, store the y value. 
        %if more than one value exists, take the mean  y value... 

        midline = zeros(size(mp2));
        hmid = round(size(midline,1)/2);
        for x = 1:size(mp2,2)
            temp = low(:,x);
            if sum(temp)== 0
                midline(hmid,x) = 1;

            elseif sum(temp) ==1
                midline(:,x) = temp;

            elseif sum(temp) >1
                idx = find(temp ==1);
                idx2 = round(median(idx));
                midline(idx2,x) = 1;
            end
        end

        %need to develop some sort of weighting. Find the average middle index
        %where the brain is (to avoid biasing from ends. 
        col_sum = sum(mp2,1);
        fst_clm_idx = find(col_sum>0,1,'first'); %number corresponds to window size
        lst_clm_idx = find(col_sum>0,1,'last'); %don't really need yet
        idx = find(midline(:,fst_clm_idx:lst_clm_idx) ==1);
        [I1,~] = ind2sub(size(mp2),idx);
        med = round(median(I1)); %median looks better than mean, less biased by outliers

        %now in each column find idx, find lowest value closest to median?
        %at a minimum, the diff should be positive
        mid2 = zeros(size(mp2));
        for x = fst_clm_idx:lst_clm_idx %just where the brain is
            if mp2(med,x) == 0;
               mid2(med,x) = 1; 

            elseif mp2(med,x) ~= 0;
                temp = mp2(:,x);
                der = diff(temp);
                idx = find(der >= 0);
                tmp = abs(idx-med);
                [~, newidx] = min(tmp); %index of closest value to estimated center
                closest = idx(newidx);
                mid2(closest,x) = 1; 
            end
        end

        %now make sure it is connected (smooth line)
        mid3 = mid2;

        for x = fst_clm_idx:lst_clm_idx-1 % assumes first one is ok
        %for each column,
            temp1 = mid3(:,x); %ref idx
            idx1 = find(temp1 ==1);
            temp2 = mid3(:,x+1); %check
            idx2 = find(temp2 ==1);
        %is it connected on either side?
            if abs(idx1-idx2) >1 
                delt = idx1-idx2; 
                if delt>0
                    idx3 = idx1 -1;
                elseif delt<0
                    idx3 = idx1 +1;
                end
                mid3(idx2,x+1) = 0;
                mid3(idx3,x+1) = 1;
            end  
        end

        %recenter first
        mid4 = zeros(img_size(1),img_size(2));
        mid4(st_idx:end_idx,:) = mid3;

        bsep = zeros(size(mid4));
        for x = 1:size(mid4,2)
            if x< fst_clm_idx %for pixels before brain
                mididx =med + st_idx;
                bsep(1:mididx,x) = 1; 
                bsep((mididx+1:end),x) = 2;      
            elseif (x>=fst_clm_idx && x<=lst_clm_idx) %for pixels in brain portion
                temp1 = mid4(:,x); %ref idx
                mididx = (find(temp1 ==1));
                bsep(1:mididx,x) = 1; 
                bsep((mididx+1:end),x) = 2;
            elseif x > lst_clm_idx
                mididx =med + st_idx;
                bsep(1:mididx,x) = 1; 
                bsep((mididx+1:end),x) = 2;
            end
        end
       %store values that are outside the area in question
        stored_vals(1:st_idx,:,i) =bsep(1:st_idx,:);
        stored_vals(end_idx:end,:,i) =bsep(end_idx,end,:);
        
    else %if there are not enough data points
        tempseg = squeeze(finalseg(:,:,i));
        %anything above mid
        bsep = zeros(img_size(1),img_size(2));
        bsep(1:(half_size(1)-1),:) = 1; %tempseg(1:(half_size(1)-1),:);
        %anything non zero value below mid way
        bsep(half_size(1):end,:) = 2; % 2.*tempseg(half_size(1):end,:);
        stored_vals(:,:,i) = bsep;
    end
    half_brain(:,:,i) = bsep;
end

% figure;
% imagesc(tempseg)
% figure;
% imagesc(squeeze(half_brain(:,:,136)))



IblurX1 = split_and_smooth(half_brain,[1 20 1]);
% 
% figure;
% title('check result')
% imshow3Dfull(IblurX1.*T1seg);

%in the y-direction, lets sample neighbouring pixels to ensure we have
%lowest intensity value chosen

half_brain2 = zeros(img_size);
sig2 = imgaussfilt3(input_image,2); %get rid of some of the noise


for y = 1:img_size(2)
    seg_slice = squeeze(IblurX1(:,y,:));
    img_slice = squeeze(sig2(:,y,:));

    %for each column, of segmentation
    new_seg = zeros(size(seg_slice));
    for x = 1:img_size(1)
        tempseg = seg_slice(:,x); %ref idx
        tempimg = img_slice(:,x); %ref idx
    %find the index of the change point
        idxch = diff(tempseg);
        idx1 = find(idxch >0 );
    %sample the voxel above and below in image, take a window +/- 3
        wind = zeros(size(tempseg));
        wind((idx1-3):(idx1+3)) = tempimg((idx1-3):(idx1+3));
    %take the index of the lowest intensity in the image
        if max(wind) >0 %make sure there are values
            idx2 = find(wind == min(wind(wind ~=0)));
        elseif max(wind) == 0 
            idx2 = idx1;
        end
        
        if size(idx2) >1
            idx2 = round(median(idx2));
        end
    %create new seg %could loop and keep old one, but i dont think its
    %necessary
        tempseg_n = tempseg;
        tempseg_n(1:idx2) = 1;
        tempseg_n(idx2+1:end) = 2;
    %fill in image
    new_seg(:,x) = tempseg_n;
    end
half_brain2(:,y,:) = new_seg;
end

% figure;
% imagesc(tempimg)

%lets try the smoothing in the y- direction. 

%quick and dirty directional gaussian smoothing. then threshold
IblurX2 = split_and_smooth(half_brain2,[1 10 5]);

%one more run through neighbouring pixels, with smaller window.

half_brain3 = zeros(img_size);
%%%%% 
for y = 1:img_size(2) %for each  coronal slice
    seg_slice = squeeze(IblurX2(:,y,:));
    img_slice = squeeze(sig2(:,y,:));

    %for each column running left-right of brain, of segmentation
    new_seg = zeros(size(seg_slice));
    for x = 1:img_size(1)
        tempseg = seg_slice(:,x); %ref idx
        tempimg = img_slice(:,x); %ref idx
    %find the index of the change point
        idxch = diff(tempseg);
        idx1 = find(idxch >0 );
    %sample the voxel above and below in image, take a window +/- 2
        wind = zeros(size(tempseg));
        wind((idx1-2):(idx1+2)) = tempimg((idx1-2):(idx1+2));
    %take the index of the lowest intensity in the image
        if max(wind) >0 %make sure there are values
            idx2 = find(wind == min(wind(wind ~=0)));
        elseif max(wind) == 0 
            idx2 = idx1;
        end
        
        if size(idx2) >1
            idx2 = round(median(idx2));
        end
    %create new seg 
        tempseg_n = tempseg;
        tempseg_n(1:idx2) = 1;
        tempseg_n(idx2+1:end) = 2;
    %fill in image
    new_seg(:,x) = tempseg_n;
    end
half_brain3(:,y,:) = new_seg;
end

IblurX3 = split_and_smooth(half_brain3,[1 8 8]); % 1 5 5 

%reset all values om IblurX3 with any in stored val
%have issues with one side and the bluring, so this should fix
for i = 1:size(finalseg,3) %for each slice
    if (sum(sum(mp(1:10,:,i))) >50) && (sum(sum(mp(11:21,:,i))) >50)
        IblurX3(:,:,i) = IblurX3(:,:,i);
    else
        IblurX3(:,:,i) = stored_vals(:,:,i);
    end 
end


%%%% New iteration
half_brain4 = zeros(img_size);
% one last run with a window of size 3, and non-smoothed image, different
% orientation
for y = 1:img_size(3) %for each axial slice
    %y = 95
    seg_slice = squeeze(IblurX3(:,:,y));
    img_slice = squeeze(T1seg(:,:,y));
    %for each column running left-right of brain, of segmentation
    new_seg = zeros(size(seg_slice));
    for x = 1:img_size(2)
        %x=154
        tempseg = seg_slice(:,x); %ref idx
        tempimg = img_slice(:,x); %ref idx
    %find the index of the change point
        idxch = diff(tempseg);
        idx1 = find(idxch >0 );
    %sample the voxel above and below in image, take a window +/- 1
        wind = tempimg((idx1-1):(idx1+1));
    %take the index of the lowest intensity in the window
        if max(wind) == 0 %if no non-zero values, then don't change
            idx2 = 2;           
        elseif max(wind) > 0 
            idx2 = find(wind == min(wind));
            idx2 = floor(mean(idx2)); %control for multiple values    
        end
    %adjust segmentation for lowest value
        if idx2 == 2 
            tempseg_n = tempseg;
        elseif idx2 == 1
            tempseg_n = tempseg;
            tempseg_n(idx1) = 2;
            tempseg_n(idx1-1) = 2;
        else 
            tempseg_n = tempseg;
            tempseg_n(idx1+1) = 1;
        end
    %fill in image
    new_seg(:,x) = tempseg_n;
    end
half_brain4(:,:,y) = new_seg;
end

IblurX4 = split_and_smooth(half_brain4,[1 5 5]); % 1 5 5 

% figure;
% title('check result')
% imshow3Dfull(IblurX4.*T1seg);
% 
% figure;
% title('check result')
% imshow3Dfull(IblurX3.*T1seg);

left_mask = (IblurX4==1) .* finalseg;
right_mask = (IblurX4==2) .* finalseg;
function [ threshold ] = threshold_finder_gaussfit(input_img,option)

%option:
%  CSF threshold on normalized image is at the first inflection point = "norm_csf". 
%  CSF threshold on sharpened image  = "sharp_csf"

% input_img = MRimage;  %troubleshooting
%  input_img = img_sharp;
% input_img = norm_img;
% input_img = wm_crop;

DS_2_hist = input_img;
DS_2_hist(input_img>1.3) = 0; %remove any upper end noise voxels to get a good histogram

%calculate intensity histogram
[h, edges] = histcounts(DS_2_hist(DS_2_hist>0.05),100);

if size(edges,2) ~= size(h,2)
    edges = edges(1:end-1);
end

% figure;
% plot(edges,h);

%Extract values for peaks
[gm_peak, wm_peak] = peak_fitter_mriseg(edges, h);

if strcmp(option, 'norm_csf') % use second derivative to determine min or max,
    
    gm_std = sqrt((gm_peak(2)/2));
    threshold = gm_peak(1) - gm_std;
    
elseif strcmp(option, 'gm_peak')
    
    threshold = gm_peak(1);

elseif strcmp(option, 'sharp_csf') %find the CSF-GM boundary minimum
    idx1 = find(edges < 0.15,1,'last');
    idx2 = find(edges > 0.5,1,'first');
    est_gm_csf = min(h(idx1:idx2));
    gm_csf = edges(h == est_gm_csf);
    threshold = (1 + 2*gm_csf) /3;
    
elseif strcmp(option, 'wm_peak')
    threshold = wm_peak(1);
end


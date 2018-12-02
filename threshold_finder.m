function [ threshold ] = threshold_finder(input_img,option)
%Use this to find the intensity value to which to threshold out the CSF. 

%option:
%  CSF threshold on normalized image is at the first inflection point = "norm_csf". 
%  CSF threshold on sharpened image  = "sharp_csf"


% zero crossing code using circshift taken from
% https://www.mathworks.com/matlabcentral/answers/267222-easy-way-of-finding-zero-crossing-of-a-function 

%input_img = MRimage;  %troubleshooting
% input_img = img_sharp;
%input_img = norm_img;
%input_img = adj_img;

DS_2_hist = input_img;
DS_2_hist(input_img>1.3) = 0; %remove any upper end noise voxels to get a good histogram

if strcmp(option, 'norm_csf') %lower threshold is 0.15
    [h, edges] = histcounts(DS_2_hist(DS_2_hist>0.05),100); % edges = upper bound of bin
    xx = linspace(0.05,1.3,200);
elseif strcmp(option, 'gm_peak') %lower threshold is 0.15
    [h, edges] = histcounts(DS_2_hist(DS_2_hist>0.05),100); % edges = upper bound of bin
    xx = linspace(0.05,1.3,200);
elseif strcmp(option, 'sharp_csf')
    [h, edges] = histcounts(DS_2_hist(DS_2_hist>0.05),100);
    xx = linspace(0.05,1.3,200);
elseif strcmp(option, 'wm_peak')
    [h, edges] = histcounts(DS_2_hist(DS_2_hist>0.25),100); % edges = upper bound of bin
    xx = linspace(0.25,1.3,200);
else
        msg = 'Please specify which option for threshold. See documentation.';
        error(msg)
end

edges = edges(2:end);
pp1 = splinefit(edges',h',20,4);
yy1 = ppval(pp1,xx);

%   now take derivative to find where it goes to 0. 
d1y = diff(yy1);
d2y = diff(d1y);
roots_d1y = find(d1y(:).*circshift(d1y(:), [-1 0]) <= 0);
roots_d1y = roots_d1y + 1; % it is shifted from the differentiation

% remove root if is at the end 
roots_d1y(roots_d1y >= size(d2y,2)) = []; 
xval_roots = xx(roots_d1y);
min_max_d2y = d2y(roots_d1y);

%   find WM peak - closest peak to 1 (property of the normalized image)
wm_peak = xval_roots(min_max_d2y < 0);
wm_peak_idx = find((abs(wm_peak - 1)) == (min(abs(wm_peak - 1))));
wm_peak = wm_peak(wm_peak_idx);

%   now go back to if statements, as only the WM peak is consistent across all
%   image types. These will find specific thresholds for different
%   conditions. 

if strcmp(option, 'norm_csf') % use second derivative to determine min or max,
    
    gm_peak_cond = min_max_d2y < 0 & xval_roots < wm_peak; %second deriv neg meaning peak, and less that WM paeak for intensity
    gm_peak = max(xval_roots(gm_peak_cond));
    %gm_csf_bound_cond = min_max_d2y > 0 & xval_roots < gm_peak; %didn't
    %work if there was a second GM peak
    
    gm_csf_bound_cond = min_max_d2y > 0 & xval_roots < 0.65;
    gm_csf_bound = max(xval_roots(gm_csf_bound_cond));
    
    if gm_csf_bound < 0.35 % this loop fixes if there is only a gradual increase in the histogram to the GM peak
        gm_csf_bound = (2*gm_peak + 3*gm_csf_bound) /5;
    end
    
    %want a point 1/3 of the way up from the bottom to control for partial
    %voluming. 
    threshold = (gm_peak + 2*gm_csf_bound) /3;
    
elseif strcmp(option, 'gm_peak')
    gm_peak_cond = min_max_d2y < 0 & xval_roots < wm_peak; %second deriv neg meaning peak, and less that WM paeak for intensity
    threshold = max(xval_roots(gm_peak_cond));

elseif strcmp(option, 'sharp_csf') %find the CSF-GM boundary minimum
    gm_csf_bound_cond = min_max_d2y > 0 & xval_roots < (wm_peak - 0.2); % sometimes GM peak isnt visible in sharp.
    gm_csf_bound = min(xval_roots(gm_csf_bound_cond));
    threshold = (wm_peak + 2*gm_csf_bound) /3;
    
%     idx1 = find(d1y>0);  %old way
%     threshold = xx(idx1(1));

elseif strcmp(option, 'wm_peak')
    threshold = wm_peak;
end

% 
% figure;
% plot(edges',h','.',xx,yy1);
% legend('data','fit')
% figure;
% plot(xx(2:end),d1y);



% legacy code 

% if strcmp(option, 'norm_csf')
%   gm_peak_cond = min_max_d2y < 0 & xval_roots < wm_peak; %second deriv neg meaning peak, and less that WM paeak for intensity
%     gm_peak = max(xval_roots(gm_peak_cond));
%     %gm_csf_bound_cond = min_max_d2y > 0 & xval_roots < gm_peak; %didn't
%     %work if there was a second GM peak
%     
%     gm_csf_bound_cond = min_max_d2y > 0 & xval_roots < 0.65;
%     gm_csf_bound = max(xval_roots(gm_csf_bound_cond));
%     
%     %want a point 1/3 of the way up from the bottom to control for partial
%     %voluming. 
%     threshold = (gm_peak + 2*gm_csf_bound) /3;
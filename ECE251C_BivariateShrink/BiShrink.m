%function img_den = BiShrink(img_n, wname, window_size, num_dec)
% Bivariate Shrinkage proposed in [5]
% 
% img_n = target image
% wname = wavelet type
% window_size = int N for NxN window
% num_dec = number of decomposition levels
% 
% 
% Algorithm: 
%   * img_n to be resized to 512x512
%     Wavelet Transform Decomposition
%     Parameter Calculations
%         noise variance
%         signal variance
%     Bivariate Shrinkage
%     Inverse Transform
%

%% Testing Variables %%
img_n = imread("lena_gray_512.tif");
wname = 'db8';

% for finding local variance at scale j of wavelet coefficients in decomp
window_size = 7;
global neighborhood;
neighborhood = ones(1,window_size)/window_size;

var_gauss = 0.025;
img_n = imnoise(img_n,"gaussian",0,var_gauss);

num_dec = 6;

%% Wavelet Decomposition (L=6)
% dwt2 is a single level decomposition (aka 1iter.) = finest level
[LL1,LH1,HL1,HH1] = dwt2(img_n,'db8', 'mode', 'per'); %finest
[LL2,LH2,HL2,HH2] = dwt2(LL1,'db8', 'mode', 'per');
[LL3,LH3,HL3,HH3] = dwt2(LL2,'db8', 'mode', 'per');
[LL4,LH4,HL4,HH4] = dwt2(LL3,'db8', 'mode', 'per');
[LL5,LH5,HL5,HH5] = dwt2(LL4,'db8', 'mode', 'per');
[LL6,LH6,HL6,HH6] = dwt2(LL5,'db8', 'mode', 'per'); %course


%% Parameter Calculations (single ex.)

% Noise variance at FINEST SCALE HH1 subband, yi specifically IN HH1
global var_noise
var_noise = median(abs(HH1), "all")/0.6745;

% Marginal Variance kth wavelet coeff
% yk = y1k, y2k, kth pixel has a kth coefficient
% this kth pixel has a neighborhood NxN centered around k
% We observe the kth child wavelet coefficient y1k  which HAS a parent y2k


% for wavelet coeff subband HH1, parent subband P(S) is HH2 (j+1)
% DO FOR EACH decomp level(6) -> each detail direction (3)
var_sig = conv2(neighborhood, neighborhood, HH1.^2, 'same');
std_wcoeff = sqrt(max(var_sig-var_noise, eps));
thresh = sqrt(3)*var_noise./std_wcoeff;


%% Bivariate Shrinkage (single ex.)
% expand parent (HH2) to be same as child (HH1)
HH2 = expand(HH2);
% MAP estimator
A = sqrt(abs(HH1).^2 + abs(HH2).^2);
numer = (A-thresh) > 0;
denom = sqrt(abs(HH1).^2 + abs(HH2).^2);
w1 = (numer/denom)*HH1;

% iterate process?
% inverse transform?

%% Complete Implementation


% img_n should be a square matrix
% length(img_n) should be a dyadic number

dwtmode('per','nodisplay');
if (num_dec<=0)
    num_level = log2(length(img_n)); % full decomposition
else
    num_level = num_dec;
end
levels = cell(1,num_level);
for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
    levels{i} = zeros(2^(i-1),2^(i-1),4); % a1 h1 v1 d1 ; a2 h2 v2 d2 ; a3 h3 v3 d3 ...
end

% additional cell array for found threshold values (based on local var)
thresh = cell(1,num_level);
for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
    thresh{i} = zeros(2^(i-1),2^(i-1),4); % a1 h1 v1 d1 ; a2 h2 v2 d2 ; a3 h3 v3 d3 ...
end


% [cA_3,cH_3,cV_3,cD_3] = dwt2(img_n,wname,'mode','per'); %finest level
% [cA_2,cH_2,cV_2,cD_2] = dwt2(cA_3,wname,'mode','per');
% [cA_1,cH_1,cV_1,cD_1] = dwt2(cA_2,wname,'mode','per');  %coarse level

% decomposition

% filling in, starting from finest decomp -> num_levels decomp
for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
    if (i==log2(length(img_n)))
       [cA_t,cH_t,cV_t,cD_t] = dwt2(img_n,wname);
       levels{i}(:,:,1) = cA_t; %LL
       levels{i}(:,:,2) = cH_t; %LH
       levels{i}(:,:,3) = cV_t; %HL
       levels{i}(:,:,4) = cD_t; %HH
    else
       [cA_t,cH_t,cV_t,cD_t] = dwt2(levels{i+1}(:,:,1),wname);
       levels{i}(:,:,1) = cA_t;
       levels{i}(:,:,2) = cH_t;
       levels{i}(:,:,3) = cV_t;
       levels{i}(:,:,4) = cD_t;
    end
end
  
% calculate/fill in thresh array to correspond to decomposition
for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
  thresh{i}(:,:,2) = param(levels{i}(:,:,2));
  thresh{i}(:,:,3) = param(levels{i}(:,:,3));
  thresh{i}(:,:,4) = param(levels{i}(:,:,4));
end

% apply shrinkage algorithm using found vars array
for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1            % trying to input corresponding thresh value 
                                                                        % vv application of this is probably wrong
  levels{i}(:,:,2) = BiShrinkage(levels{i}(:,:,2), levels{i-1}(:,:,2), thresh{i}(:,:,2));
  levels{i}(:,:,3) = BiShrinkage(levels{i}(:,:,3), levels{i-1}(:,:,3), thresh{i}(:,:,3));
  levels{i}(:,:,4) = BiShrinkage(levels{i}(:,:,4), levels{i-1}(:,:,4), thresh{i}(:,:,4));
end                                                    %   ^ how to access courser scale parent?

 
% for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
%         levels{i}(:,:,2) = nlfilter(levels{i}(:,:,2),[window_size window_size],f);
%         levels{i}(:,:,3) = nlfilter(levels{i}(:,:,3),[window_size window_size],f);
%         levels{i}(:,:,4) = nlfilter(levels{i}(:,:,4),[window_size window_size],f);
% end
% 
for i=log2(length(img_n))-num_level+1+1:log2(length(img_n))
        levels{i}(:,:,1) = idwt2(levels{i-1}(:,:,1),levels{i-1}(:,:,2),levels{i-1}(:,:,3),levels{i-1}(:,:,4),wname);
end
img_den = idwt2(levels{i}(:,:,1),levels{i}(:,:,2),levels{i}(:,:,3),levels{i}(:,:,4),wname);
imshow(img_den)




%% Functions

% calculation of threshold parameter based on neighborhood local variance
function thresh = param(child) %child, ie wavelet coeff in HH1 (with parent in HH2)
    global neighborhood
    global var_noise
    var_sig = conv2(neighborhood, neighborhood, child.^2, 'same'); 
    std_wcoeff = sqrt(max(var_sig-var_noise, eps));
    thresh = sqrt(3)*var_noise./std_wcoeff;
end

% BiShrink Algorithm: MAP estimator
% uses wavelet @ scale j, parent @scale j+1, and thresh value
function w = BiShrinkage(child, parent, thresh)
    parent = expand(parent);    % expand parent (HH2) to be same as child (HH1)
    
    % MAP estimator
    A = sqrt(abs(child).^2 + abs(parent).^2);
    numer = (A-thresh) > 0;
    w = (numer/A)*child;
end

%% Drafting 

% nlfilter(A,[m n],fun)
% nlfilter(img_n, [window_size window_size], fun)

% window_size = 7;
% global neighborhood;
% neighborhood = ones(1,window_size)/window_size


function thresh = thresh_calc(block)
    global neighborhood
    var_sig = conv2(neighborhood, neighborhood, block.^2, 'same'); 
    std_wcoeff = sqrt(max(var_sig-var_noise, eps));
    thresh = sqrt(3)*var_noise./std_wcoeff;

    
end

% thresh1 = thresh_calc(HH1);
% w1 = BiShrinkage(HH1, HH2, thresh1);
% 
% 
%nlfilter applies f function to A
% thresh = zeroes(size(HH1));
% 
% for i = 1:1:size(HH1)-1
%     for j = 1:1:size(HH1)-1
%         thresh(i,j) = nlfilter(levels{i}(:,:,2),[window_size window_size],thresh_calc);
%     end
% end
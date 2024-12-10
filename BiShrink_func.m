function img_den = BiShrink_func(img_n,wname,window_size,num_dec)
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

    neighborhood = ones(1,window_size)/window_size;
    var_noise = (median(abs(levels{log2(length(img_n))}(:,:,4)), "all")/0.6745)^2; % noise variance

    % calculate/fill in thresh array to correspond to decomposition
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
        thresh{i}(:,:,2) = param(levels{i}(:,:,2),neighborhood,var_noise);
        thresh{i}(:,:,3) = param(levels{i}(:,:,3),neighborhood,var_noise);
        thresh{i}(:,:,4) = param(levels{i}(:,:,4),neighborhood,var_noise);
    end

    %apply shrinkage algorithm using found vars array
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1+1            % trying to input corresponding thresh value
        % vv application of this is probably wrong
        levels{i}(:,:,2) = BiShrinkage(levels{i}(:,:,2), levels{i-1}(:,:,2), thresh{i}(:,:,2));
        levels{i}(:,:,3) = BiShrinkage(levels{i}(:,:,3), levels{i-1}(:,:,3), thresh{i}(:,:,3));
        levels{i}(:,:,4) = BiShrinkage(levels{i}(:,:,4), levels{i-1}(:,:,4), thresh{i}(:,:,4));
    end

    for i=log2(length(img_n))-num_level+1+1:log2(length(img_n))
        levels{i}(:,:,1) = idwt2(levels{i-1}(:,:,1),levels{i-1}(:,:,2),levels{i-1}(:,:,3),levels{i-1}(:,:,4),wname);
    end
    img_den = idwt2(levels{i}(:,:,1),levels{i}(:,:,2),levels{i}(:,:,3),levels{i}(:,:,4),wname);

end



%% Functions

% calculation of threshold parameter based on neighborhood local variance
function thresh = param(child, neighborhood, var_noise) %child, ie wavelet coeff in HH1 (with parent in HH2)
    % global neighborhood
    % global var_noise
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
    A = A - thresh;
    A  = A .* (A > 0);
    w = child .* A./(A+thresh); 
    % numer = (A-thresh) > 0;
    % w = (numer.*child)./A;
end

function [y] = expand(x)
% WAVELET SOFTWARE AT POLYTECHNIC UNIVERSITY, BROOKLYN, NY
% http://taco.poly.edu/WaveletSoftware/

[N,M] = size(x);
N = N*2;
M = M*2;

y = zeros(N,M);
y(1:2:N,1:2:M) = x;
y(2:2:N,2:2:M) = x;
y(1:2:N,2:2:M) = x;
y(2:2:N,1:2:M) = x;

end

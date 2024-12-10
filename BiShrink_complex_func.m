function img_den = BiShrink_complex_func(img_n,window_size,num_dec)
    % img_n should be a square matrix
    % length(img_n) should be a dyadic number
    
    dwtmode('per','nodisplay');
    if (num_dec<=0)
        num_level = log2(length(img_n)); % full decomposition
    else
        num_level = num_dec;
    end

    % decomposition
    [a,d] = dualtree2(img_n,"Level",num_level,"FilterLength",6,"LevelOneFilter","nearsym13_19");

    
    
    % additional cell array for found threshold values (based on local var)
    thresh = cell(1,num_level);
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
        thresh{i} = zeros(2^(i-1),2^(i-1),6); % hl1 hh1 lh1 lh1 hh1 hl1 ; hl2 hh2 lh2 lh2 hh2 hl2 ; ...
    end

    

    neighborhood = ones(1,window_size)/window_size;
    var_noise = (median(abs(d{1}(:,:,2)),"all")/0.6745)^2; % noise variance

    cnt = 1;
    % calculate/fill in thresh array to correspond to decomposition
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
        thresh{i}(:,:,1) = param(abs(d{cnt}(:,:,1)),neighborhood,var_noise);
        thresh{i}(:,:,2) = param(abs(d{cnt}(:,:,2)),neighborhood,var_noise);
        thresh{i}(:,:,3) = param(abs(d{cnt}(:,:,3)),neighborhood,var_noise);
        thresh{i}(:,:,4) = param(abs(d{cnt}(:,:,4)),neighborhood,var_noise);
        thresh{i}(:,:,5) = param(abs(d{cnt}(:,:,5)),neighborhood,var_noise);
        thresh{i}(:,:,6) = param(abs(d{cnt}(:,:,6)),neighborhood,var_noise);
        cnt = cnt + 1;
    end

    cnt = 1;
    %apply shrinkage algorithm using found vars array
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1+1            % trying to input corresponding thresh value
        % vv application of this is probably wrong
        d{cnt}(:,:,1) = BiShrinkage(d{cnt}(:,:,1), d{cnt+1}(:,:,1), thresh{i}(:,:,1));
        d{cnt}(:,:,2) = BiShrinkage(d{cnt}(:,:,2), d{cnt+1}(:,:,2), thresh{i}(:,:,2));
        d{cnt}(:,:,3) = BiShrinkage(d{cnt}(:,:,3), d{cnt+1}(:,:,3), thresh{i}(:,:,3));
        d{cnt}(:,:,4) = BiShrinkage(d{cnt}(:,:,4), d{cnt+1}(:,:,4), thresh{i}(:,:,4));
        d{cnt}(:,:,5) = BiShrinkage(d{cnt}(:,:,5), d{cnt+1}(:,:,5), thresh{i}(:,:,5));
        d{cnt}(:,:,6) = BiShrinkage(d{cnt}(:,:,6), d{cnt+1}(:,:,6), thresh{i}(:,:,6));
        cnt = cnt + 1;
    end

    img_den = idualtree2(a,d,"FilterLength",6,"LevelOneFilter","nearsym13_19");


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

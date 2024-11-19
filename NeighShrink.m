function img_den = NeighShrink(img_n,wname,window_size,num_dec)
    % this is the NeighShrink proposed in [6]. We used "Method 1" in it.
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


    % [cA_3,cH_3,cV_3,cD_3] = dwt2(img_n,wname,'mode','per'); %finest level
    % [cA_2,cH_2,cV_2,cD_2] = dwt2(cA_3,wname,'mode','per');
    % [cA_1,cH_1,cV_1,cD_1] = dwt2(cA_2,wname,'mode','per');  %coarse level
    % decomposition
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
        if (i==log2(length(img_n)))
           [cA_t,cH_t,cV_t,cD_t] = dwt2(img_n,wname);
           levels{i}(:,:,1) = cA_t;
           levels{i}(:,:,2) = cH_t;
           levels{i}(:,:,3) = cV_t;
           levels{i}(:,:,4) = cD_t;
        else
           [cA_t,cH_t,cV_t,cD_t] = dwt2(levels{i+1}(:,:,1),wname);
           levels{i}(:,:,1) = cA_t;
           levels{i}(:,:,2) = cH_t;
           levels{i}(:,:,3) = cV_t;
           levels{i}(:,:,4) = cD_t;
        end
    end
    
    %estimate noise level using the robust estimate proposed by Donoho
    noiselev = median(abs(levels{log2(length(img_n))}(:,:,4)),"all")/0.6745;
    global thresh;
    thresh = sqrt(2*log(length(img_n)^2))*noiselev; %universal threshold
    global window_size_g;
    window_size_g = window_size;

    %denoising
    f = @NeighOp;
    %f = @dummy;
    for i=log2(length(img_n)):-1:log2(length(img_n))-num_level+1
        levels{i}(:,:,2) = nlfilter(levels{i}(:,:,2),[window_size window_size],f);
        levels{i}(:,:,3) = nlfilter(levels{i}(:,:,3),[window_size window_size],f);
        levels{i}(:,:,4) = nlfilter(levels{i}(:,:,4),[window_size window_size],f);
    end

    for i=log2(length(img_n))-num_level+1+1:log2(length(img_n))
        levels{i}(:,:,1) = idwt2(levels{i-1}(:,:,1),levels{i-1}(:,:,2),levels{i-1}(:,:,3),levels{i-1}(:,:,4),wname);
    end
    img_den = idwt2(levels{i}(:,:,1),levels{i}(:,:,2),levels{i}(:,:,3),levels{i}(:,:,4),wname);


end

function p = NeighOp(block) %implements the per-pixel operation of NeighShrink
    s = sumsqr(block);
    global thresh;
    global window_size_g;
    thresh_opt = thresh * 0.25;
    t = 1 - ((window_size_g^2)*(thresh_opt^2)/s);
    beta = max([t,0]);
    center = (window_size_g-1)/2 + 1;
    p = block(center,center) * beta;
end

function p = dummy(block)
    global window_size_g;
    center = (window_size_g-1)/2 + 1;
    p = block(center,center);
end
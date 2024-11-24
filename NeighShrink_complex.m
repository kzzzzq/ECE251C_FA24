function img_den = NeighShrink_complex(img_n,window_size,num_dec)
    if (num_dec<=0)
        num_level = log2(length(img_n)); % full decomposition
    else
        num_level = num_dec;
    end
    % decomposition
    [a,d] = dualtree2(img_n,"Level",num_level);
    
    %estimate noise level using the robust estimate proposed by Donoho
    %[~,~,~,cD_t] = dwt2(img_n,'db4');
    %noiselev = median(abs(cD_t),"all")/0.6745;
    noiselev = median(abs(d{1}(:,:,2)),"all")/0.6745;
    global thresh;
    thresh = sqrt(2*log(length(img_n)^2))*noiselev; %universal threshold
    global window_size_g;
    window_size_g = window_size;

    %denoising
    f = @NeighOp_complex;
    %f = @dummy;
    for i=1:length(d)
        %stages = d{i};
        %subband = zeros(length(stages),length(stages));
        for j=1:6
            %subband = stages(:,:,j);
            subband = d{i}(:,:,j);
            d{i}(:,:,j) = nlfilter(subband,[window_size window_size],f);
        end
    end
    img_den = idualtree2(a,d);

end

function p = NeighOp_complex(block) %implements the per-pixel operation of NeighShrink
    %s = sum(sum(block.*conj(block)))/(length(block)^2);    
    s = sum(sum(block.*conj(block)));
    global thresh;
    global window_size_g;
    thresh_opt = thresh * 0.2;
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
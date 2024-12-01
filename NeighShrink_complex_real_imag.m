function img_den = NeighShrink_complex_real_imag(img_n,window_size,num_dec)
    if (num_dec<=0)
        num_level = log2(length(img_n)); % full decomposition
    else
        num_level = num_dec;
    end
    % decomposition
    [a,d] = dualtree2(img_n,"Level",num_level,"FilterLength",6,"LevelOneFilter","nearsym13_19");
    
    %estimate noise level using the robust estimate proposed by Donoho
    [~,~,~,cD_t] = dwt2(img_n,'db4');
    %noiselev = median(abs(cD_t),"all")/0.6745;
    %noiselev_real = median(abs(real(d{1}(:,:,2))),"all")/0.6745;
    %noiselev_imag = median(abs(imag(d{1}(:,:,2))),"all")/0.6745;
    noiselev_real = median(abs(d{1}(:,:,2)),"all")/0.6745;
    noiselev_imag = noiselev_real;
    global thresh_real;
    thresh_real = sqrt(2*log(length(img_n)^2))*noiselev_real; %universal threshold
    global thresh_imag;
    thresh_imag = sqrt(2*log(length(img_n)^2))*noiselev_imag; %universal threshold
    global window_size_g;
    window_size_g = window_size;

    %denoising
    f = @NeighOp_complex_real;
    g = @NeighOp_complex_imag;
    %f = @dummy;
    for i=1:length(d)
        %stages = d{i};
        %subband = zeros(length(stages),length(stages));
        for j=1:6
            %subband = stages(:,:,j);
            subband_real = real(d{i}(:,:,j));
            subband_imag = imag(d{i}(:,:,j));
            d{i}(:,:,j) = nlfilter(subband_real,[window_size window_size],f)+...
                          1j*nlfilter(subband_imag,[window_size window_size],g);
        end
    end
    img_den = idualtree2(a,d,"FilterLength",6,"LevelOneFilter","nearsym13_19");

end

function p = NeighOp_complex_real(block) %implements the per-pixel operation of NeighShrink
    %s = sum(sum(block.*conj(block)))/(length(block)^2);  
    % block = abs(block);
    s = sum(sum(block.*conj(block)));
    global thresh_real;
    global window_size_g;
    thresh_opt = thresh_real * 0.15;
    % thresh_opt = thresh * 0.15;
    t = 1 - ((window_size_g^2)*(thresh_opt^2)/s);
    beta = max([t,0]);
    center = (window_size_g-1)/2 + 1;
    p = block(center,center) * beta;
end

function p = NeighOp_complex_imag(block) %implements the per-pixel operation of NeighShrink
    %s = sum(sum(block.*conj(block)))/(length(block)^2);  
    % block = abs(block);
    s = sum(sum(block.*conj(block)));
    global thresh_imag;
    global window_size_g;
    thresh_opt = thresh_imag * 0.15;
    % thresh_opt = thresh * 0.15;
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
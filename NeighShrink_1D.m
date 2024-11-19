function sig_den = NeighShrink_1D(sig_n,wname,window_size,num_dec)
    % this is the 1D version of the NeighShrink proposed in [6], and is
    % equivalent to the NeighCoeff method by Cai and Silverman (2001)
    % sig_n should be a column vector
    % length(sig_n) should be a dyadic number
    dwtmode('per','nodisplay');
    if (num_dec<=0)
        num_level = log2(length(sig_n)); % full decomposition
    else
        num_level = num_dec;
    end
    levels = cell(1,num_level);
    for i=log2(length(sig_n)):-1:log2(length(sig_n))-num_level+1
        levels{i} = zeros(2^(i-1),2); % a1 d1 ; a2 d2 ; a3 d3 ...
    end
    
    % full decomposition
    
    for i=log2(length(sig_n)):-1:log2(length(sig_n))-num_level+1
        if (i==log2(length(sig_n)))
           [cA_t,cD_t] = dwt(sig_n,wname);
           levels{i}(:,1) = cA_t;
           levels{i}(:,2) = cD_t;
        else
           [cA_t,cD_t] = dwt(levels{i+1}(:,1),wname);
           levels{i}(:,1) = cA_t;
           levels{i}(:,2) = cD_t;
        end
    end
    %sig_den = levels;
    
    %estimate noise level using the robust estimate proposed by Donoho
    noiselev = median(abs(levels{log2(length(sig_n))}(:,2)),"all")/0.6745;
    global thresh;
    thresh = sqrt(2*log(length(sig_n)))*noiselev; %universal threshold
    global window_size_g;
    window_size_g = window_size;

    %denoising
    f = @Neigh1DOp;
    for i=log2(length(sig_n)):-1:log2(length(sig_n))-num_level+1
        levels{i}(:,2) = nlfilter(levels{i}(:,2),[window_size 1],f);
    end

    %inverse transform
    for i=log2(length(sig_n))-num_level+1+1:log2(length(sig_n))
        levels{i}(:,1) = idwt(levels{i-1}(:,1),levels{i-1}(:,2),wname);
    end
    sig_den = idwt(levels{i}(:,1),levels{i}(:,2),wname);
end

function p = Neigh1DOp(block)
    s = sumsqr(block);
    global thresh;
    global window_size_g;
    thresh_opt = thresh * 0.3;
    t = 1 - ((window_size_g)*(thresh_opt^2)/s);
    beta = max([t,0]);
    center = (window_size_g-1)/2 + 1;
    p = block(center) * beta;

end
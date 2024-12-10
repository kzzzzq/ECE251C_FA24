function [PSNR,SNR,ssim_den,img_den] = den_compare(img, noise_var,suppress_BM3D)
    % this function adds Gaussian noise with 0 mean and variance noise_var
    % to img and apply various denoising methods on it. The intensity of 
    % pixel in img should between 0 and 1. 
    % current methods: 
    % SUREShrink: db4, soft threshold
    % NeighShrink: db4
    % standard BM3D
    % NeighShrink using CWT that only shrinks magnitude
    % NeighShrink using CWT that shrinks both real and imag part. 
    % NeighShrink followed by standard BM3D
    % NeighShrink followed by the second stage of BM3D
    
    % Possible noise types to be generated 'gw', 'g1', 'g2', 'g3', 'g4', 'g1w',
    % 'g2w', 'g3w', 'g4w'.
    noise_type =  'gw';
    %noise_var = 0.02; % Noise variance
    seed = 0; % seed for pseudorandom noise realization
    [noise, PSD, kernel] = getExperimentNoise(noise_type, noise_var, seed, size(img));
    img_n = img + noise;

    subplot(3,2,1)
    imagesc(img);
    colormap gray
    title("original image")
    
    subplot(3,2,2)
    imagesc(img_n);
    colormap gray
    title("noisy image")

    %beforeSNR = ...
    %    20*log10(norm((img(:)))/norm((img(:))-(img_n(:)))); %SNR before denoising
    [beforePSNR,beforeSNR] = psnr(img_n,img);
    ssim_noise = ssim(img_n,img);
    
    img_den_sure = wdenoise2(img_n,3,'Wavelet','db4','DenoisingMethod','Sure','ThresholdRule','Soft','NoiseEstimate','LevelIndependent','NoiseDirection',["v","h","d"]);
    ssim_sure = ssim(img_den_sure,img); 
    
    %denoising using SUREShrink
    subplot(3,2,3)
    imagesc(img_den_sure)
    colormap gray
    title("denoised image using SUREShrink")
    %afterSNR_SURE = ...
    %   20*log10(norm(im2double(img(:)))/norm((img(:))-(img_den_sure(:)))); %SNR after SUREShrink
    
    [afterPSNR1,afterSNR1] = psnr(img_den_sure,img);
    
    
    
    img_den_nei = NeighShrink(img_n,'db4',3,4);%denoising using NeighShrink

    subplot(3,2,4)
    imagesc(img_den_nei)
    colormap gray
    title("denoised image using NeighShrink")
    %afterSNR_Neigh = ...
    %    20*log10(norm(im2double(img(:)))/norm((img(:))-(img_den_nei(:))));
    [afterPSNR2,afterSNR2] = psnr(img_den_nei,img);
    ssim_nei = ssim(img_den_nei,img);
    

    if (suppress_BM3D)
        
        afterPSNR3 = 0;
        afterSNR3 = 0;
        ssim_BM3D = 1;
    else
        subplot(3,2,5)
        % img_den_BM3D_basic = BM3D_matlab_o(img_n,sqrt(noise_var)*255); %% BM3D image denoising,basic estimation, thresholding stage.
        % img_den_BM3D = BM3D_matlab_wiener_o(img_n,img_den_BM3D_basic,sqrt(noise_var)*255);
        img_den_BM3D = BM3D(img_n, sqrt(noise_var));

        [afterPSNR3,afterSNR3] = psnr(img_den_BM3D,img);
        ssim_BM3D = ssim(img_den_BM3D,img);
        imagesc(img_den_BM3D)
        colormap gray
        title("denoised image using standard BM3D")
    end
    subplot(3,2,6)
    img_den_nei_complex = NeighShrink_complex(img_n,3,4);
    [afterPSNR4,afterSNR4] = psnr(img_den_nei_complex,img);
    ssim_nei_complex = ssim(img_den_nei_complex,img);
    imagesc(img_den_nei_complex)
    colormap gray
    title("denoised image using DTCWT NeighShrink")

    if (suppress_BM3D)
        afterPSNR5 = 0;
        afterSNR5 = 0;
        ssim_BM3D_complex = 1;
    else
        % afterPSNR5 = 0;
        % afterSNR5 = 0;
        % ssim_BM3D_complex = 1;
        figure
        subplot(3,2,1)
        %img_den_BM3D_complex = BM3D_matlab_wiener_o(img_n,img_den_nei_complex,sqrt(noise_var)*255);
        img_den_BM3D_complex = BM3D(img_n, sqrt(noise_var),'stage_arg',img_den_nei_complex);
        [afterPSNR5,afterSNR5] = psnr(img_den_BM3D_complex,img);
        ssim_BM3D_complex = ssim(img_den_BM3D_complex,img);
        imagesc(img_den_BM3D_complex)
        colormap gray
        title("denoised image using DTCWT NeighShrink + BM3D second stage")
    end


    if (suppress_BM3D)
        afterPSNR6 = 0;
        afterSNR6 = 0;
        ssim_BM3D_complex_3 = 1;
    else
        
        subplot(3,2,2)
        noise_var_c = (1/(size(img,1)^2)*sum(sum(img.*img)));
        noise_var_c = noise_var_c / (10^(afterSNR4/10));
        %img_den_BM3D_complex_basic = BM3D_matlab_o(img_den_nei_complex,sqrt(noise_var_c)*255);
        %img_den_BM3D_complex_3 = BM3D_matlab_wiener_o(img_n,img_den_BM3D_complex_basic,sqrt(noise_var)*255);
        img_den_BM3D_complex_3 = BM3D(img_den_nei_complex,sqrt(noise_var_c));
        [afterPSNR6,afterSNR6] = psnr(img_den_BM3D_complex_3,img);
        ssim_BM3D_complex_3 = ssim(img_den_BM3D_complex_3,img);
        imagesc(img_den_BM3D_complex_3)
        colormap gray
        title("denoised image using DTCWT NeighShrink + BM3D")
    end

    subplot(3,2,3)
    img_den_nei_complex_real_imag = NeighShrink_complex_real_imag(img_n,3,4);
    [afterPSNR7,afterSNR7] = psnr(img_den_nei_complex_real_imag,img);
    ssim_nei_complex_real_imag = ssim(img_den_nei_complex_real_imag,img);
    imagesc(img_den_nei_complex_real_imag)
    colormap gray
    title("denoised image using DTCWT NeighShrink (seperate real and imag parts)")

   
    subplot(3,2,4)
    img_den_bi = BiShrink_func(img_n,'db4',5,5);
    [afterPSNR8,afterSNR8] = psnr(img_den_bi,img);
    ssim_bi = ssim(img_den_bi,img);
    imagesc(img_den_bi)
    colormap gray
    title("bishrink")

    subplot(3,2,5)
    img_den_bi_com = BiShrink_complex_func(img_n,5,5);
    [afterPSNR9,afterSNR9] = psnr(img_den_bi_com,img);
    ssim_bi_com = ssim(img_den_bi_com,img);
    imagesc(img_den_bi_com)
    colormap gray
    title("ocmplex bishrink")
    
      
    disp(["noise variance:",num2str(sqrt(noise_var)*255)])
    disp(["noisy image PSNR/SSIM:",num2str(beforePSNR)," / ",num2str(ssim_noise)])
    disp(["PSNR/SSIM after SURE:",num2str(afterPSNR1)," / ",num2str(ssim_sure)])
    disp(["PSNR/SSIM after NeighShrink:",num2str(afterPSNR2)," / ",num2str(ssim_nei)])
    if (suppress_BM3D)
        
    else
        disp(["PSNR/SSIM after standard BM3D:",num2str(afterPSNR3)," / ",num2str(ssim_BM3D)])
        disp(["PSNR/SSIM after CWT Neighshrink + BM3D second stage:",num2str(afterPSNR5)," / ",num2str(ssim_BM3D_complex)])
        disp(["PSNR/SSIM after CWT Neighshrink + BM3D:",num2str(afterPSNR6)," / ",num2str(ssim_BM3D_complex_3)])
    end
    
    disp(["PSNR/SSIM after CWT NeighShrink:",num2str(afterPSNR4)," / ",num2str(ssim_nei_complex)])
    disp(["PSNR/SSIM after CWT NeighShrink (seperate real and imag parts):",num2str(afterPSNR7)," / ",num2str(ssim_nei_complex_real_imag)])
    disp(["PSNR/SSIM after BiShrink:",num2str(afterPSNR8)," / ",num2str(ssim_bi)])
    disp(["PSNR/SSIM after CWT BiShrink:",num2str(afterPSNR9)," / ",num2str(ssim_bi_com)])

    % disp(["PSNR/SSIM after CWT Neighshrink + BM3D:",num2str(afterPSNR6)," / ",num2str(ssim_BM3D_complex_3)])

    SNR = [beforeSNR,afterSNR1,afterSNR2,afterSNR3,afterSNR4,afterSNR5,afterSNR6,afterSNR7,afterSNR8,afterSNR9];
    PSNR = [beforePSNR,afterPSNR1,afterPSNR2,afterPSNR3,afterPSNR4,afterPSNR5,afterPSNR6,afterPSNR7,afterPSNR8,afterPSNR9];
    ssim_den = [ssim_noise,ssim_sure,ssim_nei,ssim_BM3D,ssim_nei_complex,ssim_BM3D_complex,ssim_BM3D_complex_3,ssim_nei_complex_real_imag,ssim_bi,ssim_bi_com];
    img_den(:,:,1) = img;
    img_den(:,:,2) = img_n;
    img_den(:,:,3) = img_den_sure;
    img_den(:,:,4) = img_den_nei;
    img_den(:,:,5) = img_den_BM3D;
    img_den(:,:,6) = img_den_nei_complex;
    img_den(:,:,7) = img_den_BM3D_complex;
    img_den(:,:,8) = img_den_BM3D_complex_3;
    img_den(:,:,9) = img_den_nei_complex_real_imag;
    img_den(:,:,10) = img_den_bi;
    img_den(:,:,11) = img_den_bi_com;
end
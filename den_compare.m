function [PSNR,SNR,ssim_den] = den_compare(img, noise_var)
    % this function adds Gaussian noise with 0 mean and variance noise_var
    % to img and apply various denoising methods on it. The intensity of 
    % pixel in img should between 0 and 1. 
    img_n = imnoise(img,"gaussian",0,noise_var); %adding noise
    
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
    

    subplot(3,2,5)
    img_den_BM3D_basic = BM3D_matlab_o(img_n,sqrt(noise_var)*255); %% BM3D image denoising,basic estimation, thresholding stage.
    img_den_BM3D = BM3D_matlab_wiener_o(img_n,img_den_BM3D_basic,sqrt(noise_var)*255);
    [afterPSNR3,afterSNR3] = psnr(img_den_BM3D,img);
    ssim_BM3D = ssim(img_den_BM3D,img);
    imagesc(img_den_BM3D)
    colormap gray
    title("denoised image using standard BM3D")

    subplot(3,2,6)
    img_den_nei_complex = NeighShrink_complex(img_n,3,4);
    [afterPSNR4,afterSNR4] = psnr(img_den_nei_complex,img);
    ssim_nei_complex = ssim(img_den_nei_complex,img);
    imagesc(img_den_nei_complex)
    colormap gray
    title("denoised image using DTCWT NeighShrink")

    figure
    subplot(3,2,1)
    img_den_BM3D_complex = BM3D_matlab_wiener_o(img_n,img_den_nei_complex,sqrt(noise_var)*255);
    [afterPSNR5,afterSNR5] = psnr(img_den_BM3D_complex,img);
    ssim_BM3D_complex = ssim(img_den_BM3D_complex,img);
    imagesc(img_den_BM3D_complex)
    colormap gray
    title("denoised image using DTCWT NeighShrink followed by BM3D")
    
    disp(["noisy image PSNR/SSIM:",num2str(beforePSNR)," / ",num2str(ssim_noise)])
    disp(["PSNR/SSIM after SURE:",num2str(afterPSNR1)," / ",num2str(ssim_sure)])
    disp(["PSNR/SSIM after NeighShrink:",num2str(afterPSNR2)," / ",num2str(ssim_nei)])
    disp(["PSNR/SSIM after standard BM3D:",num2str(afterPSNR3)," / ",num2str(ssim_BM3D)])
    disp(["PSNR/SSIM after CWT NeighShrink:",num2str(afterPSNR4)," / ",num2str(ssim_nei_complex)])
    disp(["PSNR/SSIM after CWT Neighshrink + BM3D second stage:",num2str(afterPSNR5)," / ",num2str(ssim_BM3D_complex)])


    SNR = [beforeSNR,afterSNR1,afterSNR2,afterSNR3,afterSNR4,afterSNR5];
    PSNR = [beforePSNR,afterPSNR1,afterPSNR2,afterPSNR3,afterPSNR4,afterPSNR5];
    ssim_den = [ssim_noise,ssim_sure,ssim_nei,ssim_BM3D,ssim_nei_complex,ssim_BM3D_complex];
end
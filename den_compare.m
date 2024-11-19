function [PSNR,SNR] = den_compare(img, noise_var)
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
    img_den_sure = wdenoise2(img_n,3,'Wavelet','db4','DenoisingMethod','Sure','ThresholdRule','Soft','NoiseEstimate','LevelIndependent','NoiseDirection',["v","h","d"]);
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

    subplot(3,2,5)
    img_den_BM3D_basic = BM3D_matlab(img_n,sqrt(noise_var)*255); %% BM3D image denoising,basic estimation, thresholding stage.
    img_den_BM3D = BM3D_matlab_wiener(img_n,img_den_BM3D_basic,sqrt(noise_var)*255);
    [afterPSNR3,afterSNR3] = psnr(img_den_BM3D,img);
    imagesc(img_den_BM3D)
    colormap gray
    title("denoised image using standard BM3D")

    SNR = [beforeSNR,afterSNR1,afterSNR2,afterSNR3];
    PSNR = [beforePSNR,afterPSNR1,afterPSNR2,afterPSNR3];
end
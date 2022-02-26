function output_sentil()

clear all;
close all;

addpath ./flow-code-matlab;
addpath(genpath('./csh'));
addpath ./siftDemoV4;
addpath ./utils;

all_img_path = 'C:\study\DATA\MPI-Sintel-complete\training\clean\';
imgpath = dir(all_img_path); imgpath = imgpath(3:end);
all_data_path = 'C:\study\DATA\MPI-Sintel-complete\training\flow\';
datapath = dir(all_data_path); datapath = datapath(3:end);

csh_wid = 16;
homo_num = 32;
scale = 1;

for fold_id = 8:length(imgpath)
    imglist = dir([all_img_path imgpath(fold_id).name '\*.png']);
    for ind = 1 : length(imglist)-1


        %% ground truth
        flo = readFlowFile( [all_data_path datapath(fold_id).name '\' imglist(ind).name(1:end-4) '.flo']);
        testimggt = flowToColor(flo); %figure,imshow(testimggt)
        maxrad = scale_flow(flo,50);

        filename = [all_img_path imgpath(fold_id).name '\' imglist(ind).name(1:end-4) '_knn70.jpg'];
        testimg = imread(filename);
        
        filename = [all_img_path imgpath(fold_id).name '\' imglist(ind).name(1:end-4) '_nl.jpg'];
        testimg1 = imread(filename);
        
        filename = [all_img_path imgpath(fold_id).name '\' imglist(ind).name(1:end-4) '_mdpof.flo'];
        uv2 = readFlowFile(filename);
        maxrad2 = sqrt(uv2(:,:,1).^2+uv2(:,:,2).^2);        maxrad_final = max(maxrad,maxrad2);
        tuv = uv2 ./ repmat(maxrad_final,[1 1 2]);      testimg2 = flowToColor(tuv); %subplot(223),imshow(testimg2)
        %subplot(224),imshow(testimggt)

        [m,n,c]=size(flo);
        newimg = zeros(2*m+15,2*n+15,3);
        newimg = uint8(newimg);
        newimg(6:m+5,6:n+5,:) = testimggt;
        newimg(6:m+5,n+11:2*n+10,:) = testimg;
        newimg(m+11:2*m+10,6:n+5,:) = testimg1;
        newimg(m+11:2*m+10,n+11:2*n+10,:) = testimg2;

        filename=[all_img_path imgpath(fold_id).name '\' imglist(ind).name(1:end-4) '_compare.jpg'];
        imwrite(newimg,filename);


        %% final refinement by optical flow

        close all;
    end
end
;
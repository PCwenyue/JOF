function test_flow()

clear all;
close all;

addpath ./flow-code-matlab;
addpath(genpath('./csh'));
addpath ./siftDemoV4;
addpath ./utils;
addpath ./others;

a{1}='Dimetrodon';
a{2}='Grove2';
a{3}='Grove3';
a{4}='Hydrangea';
a{5}='RubberWhale';
a{6}='Urban2';
a{7}='Urban3';
a{8}='Venus';
a{9}='Beanbags';
a{10}='DogDance';
a{11}='MiniCooper';
a{12}='Walking';
imgpath = 'E:\葛利跃\调参用\陈卓远2015 - 副本\Knn_stat\other-data\';
datapath = 'F:\葛利跃\论文程序调试好的葛利跃\classic-hs\data\other-gt-flow\';

csh_wid = 16;

for ind = 8
    %I0 = imread('tmp2.jpg');I1=imread('tmp3.jpg');
    I0 = imread([imgpath a{ind} '\frame10.png']); 
    I1 = imread([imgpath a{ind} '\frame11.png']);
    
    %% ground truth
    flo=readFlowFile( [datapath a{ind} '\flow10.flo']);
    testimggt = flowToColor(flo);
    filename=[imgpath a{ind} '\' a{ind} 'gt.png'];
    imwrite(testimggt,filename);
    maxrad = scale_flow(flo,50);
    
%      [histmat_gt top_uv1] = plot_2d2(flo,20,1); % 原
    [histmat_gt top_uv1] = plot_2d2(flo,1);
     
     %% direct CSH
     [uv1,flag1] = CSH_nn_flow_base(I0,I1,[],csh_wid,1);
%      [histmat_csh top_uv2] =
%      plot_2d2(uv1(1:end-csh_wid+1,1:end-csh_wid+1,:),20,0); % 原
  [histmat_csh top_uv2] = plot_2d2(uv1(1:end-csh_wid+1,1:end-csh_wid+1,:),20);
     
     flag2 = csh_valid(uv1,top_uv2);
     flag = flag1 .* flag2; figure,imshow(flag)

%      maxrad2 = sqrt(uv1(:,:,1).^2+uv1(:,:,2).^2);
%      maxrad_final = max(maxrad,maxrad2);
%      uv1 = uv1 ./ repmat(maxrad_final,[1 1 2]);
%      testimg1 = flowToColor(uv1); figure,imshow(testimg1)
    
    % collect uv2 into grids
    %uv_core=collect_uv(top_uv2);
    
     figure(4),  plot(top_uv1(1,:),top_uv1(2,:),'bo', 'MarkerSize',13)
     hold;    plot(top_uv2(1,:),top_uv2(2,:),'rx', 'MarkerSize',13)

    %% top homography computation by SIFT
%     top_homo1 = sift_match(I0,I1);
%     if length(top_homo1) > 30
%         top_homo1 = homo_pick(top_homo1);
%     end
    %demo_flow(top_homo1,m,n);
    
    top_homo2 = dense_homo(uv1(1:end-csh_wid,1:end-csh_wid,:),csh_wid/2-.5,flag1(1:end-csh_wid,1:end-csh_wid));
    
    close all;
    [uv,uv1,uv2] = knn_flow(I0,I1,top_uv2,top_homo2);
    
    maxrad2 = sqrt(uv(:,:,1).^2+uv(:,:,2).^2);        maxrad_final = max(maxrad,maxrad2);
    tuv = uv ./ repmat(maxrad_final,[1 1 2]);        testimg = flowToColor(tuv); figure, subplot(221),imshow(testimg)

    maxrad2 = sqrt(uv1(:,:,1).^2+uv1(:,:,2).^2);        maxrad_final = max(maxrad,maxrad2);
    tuv = uv1 ./ repmat(maxrad_final,[1 1 2]);        testimg1 = flowToColor(tuv); subplot(222),imshow(testimg1)

    maxrad2 = sqrt(uv2(:,:,1).^2+uv2(:,:,2).^2);        maxrad_final = max(maxrad,maxrad2);
    tuv = uv2 ./ repmat(maxrad_final,[1 1 2]);      testimg2 = flowToColor(tuv); subplot(223),imshow(testimg2)
    subplot(224),imshow(testimggt)
    
    imwrite(testimg,[imgpath a{ind} '\' a{ind} 'knn_fusion2.png']);
    imwrite(testimg1,[imgpath a{ind} '\' a{ind} 'knn_offset2.png']);
    imwrite(testimg2,[imgpath a{ind} '\' a{ind} 'knn_rigid2.png']);
    
    [m,n,c]=size(I0);
    newimg = zeros(2*m+15,2*n+15,3);
    newimg = uint8(newimg);
    newimg(6:m+5,6:n+5,:) = testimggt;
    newimg(6:m+5,n+11:2*n+10,:) = testimg;
    newimg(m+11:2*m+10,6:n+5,:) = testimg1;
    newimg(m+11:2*m+10,n+11:2*n+10,:) = testimg2;
    
    filename=[imgpath a{ind} '\knn_compo.png'];
    imwrite(newimg,filename);
    
    %filename=[imgpath a{ind} '\knn_uv.mat'];
    %save(filename,'uv2');
    
    %% final refinement by optical flow
   
    close all;
end
;

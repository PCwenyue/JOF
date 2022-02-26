function test_stereo2()

clear all;
close all;

addpath ./flow-code-matlab;
% addpath ./flow_code_v2;
% addpath (genpath('./flow_code_v2/utils'))
%addpath(genpath('./csh'));

imgpath = 'D:\zhuchen\data\fog_stereo\';
imglist = dir(imgpath);

for ind = 3 : length(imglist)
    I0 = imread([imgpath imglist(ind).name '\view1.png']); 
    I1 = imread([imgpath imglist(ind).name '\view5.png']);
        
    % deqing's method
    uv = estimate_flow_interface(I0, I1, 'classic+nl');
    u = -uv(:,:,1); u(u<0)=0;
    filename=[imgpath imglist(ind).name '\dis_deqing.png'];
    imwrite(uint8(u*3),filename);
        
    save([imgpath imglist(ind).name '\result2.mat'],'u');

end


function test_csh()

clear all;
close all;

addpath ./flow-code-matlab;
addpath(genpath('./csh'));

datapath = 'D:/zhuchen/data/CSH_dataset/';
datalist = dir(datapath);

ope.u_max = 100; ope.u_min = -100;
ope.v_max = 100; ope.v_min = -100;

for data_id = 3 : length(datalist)
    imglist = dir([datapath datalist(data_id).name '/*.jpg']);
    for img_id = 1:2:length(imglist)-1
        I0 = imread([datapath datalist(data_id).name '/' imglist(img_id).name]); 
        I1 = imread([datapath datalist(data_id).name '/' imglist(img_id+1).name]);
        
        [m,n,c] = size(I0);
        
        while max(m,n) > 600
            I0 = imresize(I0,.5);
            I1 = imresize(I1,.5);
            [m,n,c] = size(I0);
        end

        uv1 = CSH_nn_flow_base(I0,I1,100,8,1);
        testimg1 = flowToColor(uv1); %imshow(testimg)
        filename = [datapath datalist(data_id).name '/' imglist(img_id).name(1:end-4) 'flow1.png'];
        imwrite(testimg1,filename);

        uv2 = knn_flow(I0,I1,[],100,ope); % after multi-scale fusion-based optimization
        testimg2 = flowToColor(uv2);
        filename = [datapath datalist(data_id).name '/' imglist(img_id).name(1:end-4) 'flow2.png'];
        imwrite(testimg2,filename);
%         figure,subplot(121),imshow(testimg1)
%         subplot(122),imshow(testimg2)
        
        filename = [datapath datalist(data_id).name '/' imglist(img_id).name(1:end-4) 'flow.mat'];
        save(filename,'uv1','uv2');
    end
end

end


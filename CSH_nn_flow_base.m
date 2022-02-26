function [uv,valid_flag] = CSH_nn_flow_base(im1,im2,maxflow,patch_wid,Knn_count)
%addpath(genpath('./csh'));

[m,n,c] = size(im1);
if ~exist('maxflow','var') || isempty(maxflow)
    maxflow = max(m,n);
end

if ~exist('patch_wid','var')
    patch_wid = 8; %patch: 8 x 8
end
if ~exist('Knn_count','var')
    Knn_count = 1; % 原始为1
end

iterations = 7;display = 0; % 原始迭代7次
valid_flag = ones(m-patch_wid+1,n-patch_wid+1);

%% CSH NNF
[CSH_ann, CSH_bnn] = CSH_nn(im1,im2,... % basic params - can handle only these (rest can be default)
    patch_wid,iterations,Knn_count,1);
CSH_ann = double(CSH_ann); CSH_ann = CSH_ann(1:end-patch_wid+1,1:end-patch_wid+1,:);
CSH_bnn = double(CSH_bnn); CSH_bnn = CSH_bnn(1:end-patch_wid+1,1:end-patch_wid+1,:);
sz = size(CSH_ann);

%% move to center
[X,Y] = meshgrid(1:n-patch_wid+1,1:m-patch_wid+1);
u = CSH_ann(:,:,1) - X; u = moveToCenter(u,patch_wid);
v = CSH_ann(:,:,2) - Y; v = moveToCenter(v,patch_wid);
uv = cat(3,u,v);
valid_flag = moveToCenter(valid_flag,patch_wid);
%valid_flag(1:1+patch_wid/2:end-patch_wid/2+1,1:1+patch_wid/2:end-patch_wid/2+1) = 1;

%% bi-direction check
CSH_x1back = zeros(sz(1),sz(2),2);
for i=1:sz(1)
    for j=1:sz(2)
        tmp_x = CSH_ann(i,j,1);        tmp_y = CSH_ann(i,j,2);
        CSH_x1back(i,j,1) = CSH_bnn(tmp_y,tmp_x,1); % X back
        CSH_x1back(i,j,2) = CSH_bnn(tmp_y,tmp_x,2); % Y back
    end
end
CSH_x1back = CSH_x1back-cat(3,X,Y);
tmp_flag = (abs(CSH_x1back(:,:,1))<=1) & (abs(CSH_x1back(:,:,1))<=1);
tmp_flag = moveToCenter(tmp_flag,patch_wid);
valid_flag = tmp_flag & valid_flag;

%% maximum velocity
tmp_flag = (abs(uv(:,:,1))<maxflow) & (abs(uv(:,:,2))<maxflow);
valid_flag = tmp_flag & valid_flag;

%uv(uv>maxflow)=0;%maxflow;%uv(uv<-maxflow)=0;%-maxflow;

if display == 1
    img = flowToColor(uv);    imshow(img)
end

%bidr_flag = (CSH_x1back(:,:,1)==0).*(CSH_x1back(:,:,2)==0);

end

function Corr = moveToCenter(DCF, patch_w)
sz = size(DCF);
if length(sz) == 2
    sz = [sz 1];
end
Corr=zeros(sz(1)+patch_w-1, sz(2)+patch_w-1, sz(3));
Corr(patch_w/2+1:end-patch_w/2+1, patch_w/2+1:end-patch_w/2+1, :) = DCF(:,:,:);
end
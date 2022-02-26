function [epe ae epemap aemap]=eva_flow2(flow1,flow2,im1,im2)
%flow1: estimated flow
%flow2: ground truth

h=size(flow1,1);w=size(flow2,2);

%ae = zeros(h,w);
%epe = zeros(h,w);
hidden = zeros(h,w);

%% EPE computing
hid_id=find(flow2>h);
flow1(hid_id)=0;
flow2(hid_id)=0;
epe=(flow1-flow2).^2;
epe=epe(:,:,1)+epe(:,:,2);
epe=sqrt(epe);
epemap=epe;

%% AE computing
tmp1=1+flow1(:,:,1).*flow2(:,:,1)+flow1(:,:,2).*flow2(:,:,2);
tmp2=sqrt(1+flow1(:,:,1).^2+flow1(:,:,2).^2);
tmp3=sqrt(1+flow2(:,:,1).^2+flow2(:,:,2).^2);
ae=tmp1./tmp2./tmp3;
ae=acos(ae);
aemap=ae;

ae = sum(sum(ae))/h/w*180/pi;
epe = sum(sum(epe))/h/w;

%% show results
if exist('im1','var')
    figure,
    subplot(231),imshow(im1);
    subplot(232),imshow(im2);
    subplot(233),imshow(hidden);

%     subplot(334),quiver(flow1(1:4:end,1:4:end,1),flow1(1:4:end,1:4:end,2))
%     subplot(335),quiver(flow2(1:4:end,1:4:end,1),flow2(1:4:end,1:4:end,2))

    subplot(234),imshow(abs(ae),[])
    subplot(235),imshow(abs(epe),[])
end
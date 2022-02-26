function final_uv = qpbo_fusion(I0,I1,color_I0,uv1,uv2,this)
% use QPBO method to fuse SIFT flow and optical flow
% uv1: fused discrete flow from SIFT
% uv2: resampled continuous flow
% this.lambda: weights of spatial term
% this.

[m,n,c] = size(color_I0);
[x,y] = meshgrid(1:n,1:m);

amp_factor = 100;
threshold = 5;

if c == 3
    [weightx weighty] = weighted_derivative_ops_color(m, n, color_I0, 3, 0.01, 0);
else
    [weightx weighty] = weighted_derivative_ops_grayscale(m, n, color_I0, 3, 0.01, 0);
end
[m,n,c] = size(I0);
% error map computing for uv1
x1 = x + uv1(:,:,1);
y1 = y + uv1(:,:,2);
warpIm1 = zeros(m,n,c);
for j=1:c
    warpIm1(:,:,j) = interp2(x,y,I1(:,:,j),x1,y1,'bicubic');
end
error_map1 = I0 - warpIm1;
error_map1(isnan(error_map1))=0;
%error_map1 = inpaintnan(error_map1);
error_map1 = evaluate(this.rho_spatial_u{1}, error_map1) ;
error_map1 = mean(error_map1,3);
error_map1 = error_map1(:);

% error map computing for uv1
x2 = x + uv2(:,:,1);
y2 = y + uv2(:,:,2);
warpIm2 = zeros(m,n,c);
for j=1:c
    warpIm2(:,:,j) = interp2(x,y,I1(:,:,j),x2,y2,'bicubic');
end
error_map2 = I0 - warpIm2;
error_map2(isnan(error_map2))=0;
%error_map2 = inpaintnan(error_map2);
error_map2 = evaluate(this.rho_spatial_u{1}, error_map2) ;
error_map2 = mean(error_map2,3);
error_map2 = error_map2(:);

% QPBO fusion
UE = [evaluate(this.rho_spatial_u{1}, error_map1)     evaluate(this.rho_spatial_u{1}, error_map2)]; % 2xM unary data consistency
UE = UE';

PI = edges4connected(m,n)'; % pairwise indices, first y-direction, then x
%PE = zeros(4,size(PI,2)); % E(00,01,10,11)

%% y-direction difference
du1dy = uv1(1:end-1,:,1) - uv1(2:end,:,1);  du1dy(abs(du1dy)<threshold)=0;
dv1dy = uv1(1:end-1,:,2) - uv1(2:end,:,2);  dv1dy(abs(dv1dy)<threshold)=0;
du1dy = weighty.*du1dy(:);       dv1dy = weighty.*dv1dy(:);
E00y = evaluate(this.rho_spatial_u{1},du1dy(:)) + evaluate(this.rho_spatial_u{1},dv1dy(:)); % all flow come from uv1

du2dy = uv2(1:end-1,:,1) - uv2(2:end,:,1);  du2dy(abs(du2dy)<threshold)=0;
dv2dy = uv2(1:end-1,:,2) - uv2(2:end,:,2);  dv2dy(abs(dv2dy)<threshold)=0;
du2dy = weighty.*du2dy(:);       dv2dy = weighty.*dv2dy(:);
E11y = evaluate(this.rho_spatial_u{1},du2dy(:)) + evaluate(this.rho_spatial_u{1},dv2dy(:)); % all flow come from uv2

du12dy = uv1(1:end-1,:,1) - uv2(2:end,:,1); du12dy(abs(du12dy)<threshold)=0;
dv12dy = uv1(1:end-1,:,2) - uv2(2:end,:,2); dv12dy(abs(dv12dy)<threshold)=0;
du12dy = weighty.*du12dy(:);       dv12dy = weighty.*dv12dy(:);
E01y = evaluate(this.rho_spatial_u{1},du12dy(:)) + evaluate(this.rho_spatial_u{1},dv12dy(:)); % uv1, uv2 pairwise

du12dy = uv2(1:end-1,:,1) - uv1(2:end,:,1); du12dy(abs(du12dy)<threshold)=0;
dv12dy = uv2(1:end-1,:,2) - uv1(2:end,:,2); dv12dy(abs(dv12dy)<threshold)=0;
du12dy = weighty.*du12dy(:);       dv12dy = weighty.*dv12dy(:);
E10y = evaluate(this.rho_spatial_u{1},du12dy(:)) + evaluate(this.rho_spatial_u{1},dv12dy(:)); % uv1, uv2 pairwise

%% x-direction difference
du1dx = uv1(:,1:end-1,1) - uv1(:,2:end,1);  du1dx(abs(du1dx)<threshold)=0;
dv1dx = uv1(:,1:end-1,2) - uv1(:,2:end,2);  dv1dx(abs(dv1dx)<threshold)=0;
du1dx = weightx.*du1dx(:);       dv1dx = weightx.*dv1dx(:);
E00x = evaluate(this.rho_spatial_u{1},du1dx(:)) + evaluate(this.rho_spatial_u{1},dv1dx(:)); % all flow come from uv1

du2dx = uv2(:,1:end-1,1) - uv2(:,2:end,1);  du2dx(abs(du2dx)<threshold)=0;
dv2dx = uv2(:,1:end-1,2) - uv2(:,2:end,2);  dv2dx(abs(dv2dx)<threshold)=0;
du2dx = weightx.*du2dx(:);       dv2dx = weightx.*dv2dx(:);
E11x = evaluate(this.rho_spatial_u{1},du2dx(:)) + evaluate(this.rho_spatial_u{1},dv2dx(:)); % all flow come from uv2

du12dx = uv1(:,1:end-1,1) - uv2(:,2:end,1); du12dx(abs(du12dx)<threshold)=0;
dv12dx = uv1(:,1:end-1,2) - uv2(:,2:end,2); dv12dx(abs(dv12dx)<threshold)=0;
du12dx = weightx.*du12dx(:);       dv12dx = weightx.*dv12dx(:);
E01x = evaluate(this.rho_spatial_u{1},du12dx(:)) + evaluate(this.rho_spatial_u{1},dv12dx(:)); % uv1, uv2 pairwise

du12dx = uv2(:,1:end-1,1) - uv1(:,2:end,1); du12dx(abs(du12dx)<threshold)=0;
dv12dx = uv2(:,1:end-1,2) - uv1(:,2:end,2); dv12dx(abs(dv12dx)<threshold)=0;
du12dx = weightx.*du12dx(:);       dv12dx = weightx.*dv12dx(:);
E10x = evaluate(this.rho_spatial_u{1},du12dx(:)) + evaluate(this.rho_spatial_u{1},dv12dx(:)); % uv1, uv2 pairwise

%% PE integration
%y-direction
PEy = zeros(4,length(E00y));
PEy(1,1:length(E00y)) = E00y';  PEy(4,1:length(E11y)) = E11y';    PEy(2,1:length(E01y)) = E01y';  PEy(3,1:length(E01y)) = E10y';
%x-direction
PEx = zeros(4,length(E00x));
PEx(1,1:length(E00x)) = E00x';  PEx(4,1:length(E11x)) = E11x';    PEx(2,1:length(E01x)) = E01x';  PEx(3,1:length(E01x)) = E10x';
PE = [PEy PEx];

PE = PE * this.lambda;

% QPBO
L = vgg_qpbo(UE * amp_factor, uint32(PI), PE * amp_factor);
% eliminate some noisy labeling
L( L < 0 ) = 0;   L( L > 1 ) = 1;

% reorganize final_uv
L = reshape(L, [m n]); L = repmat(L, [1 1 2]); L = double(L);
final_uv = uv1 .* (1-L) + uv2 .* L;

end
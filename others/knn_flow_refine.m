function uv=knn_flow_refine(im1,im2,uv0,patch_wid,func)
% im1(x)=im2(x+u)
% im1, im2: images
% uv: upscaled kd-tree-flow from last pyramid
% maxflow: maximum flow 
% patch_rad: size of the patch
addpath ./Ann;
%% parameters setting
if isa(im1,'uint8')
    im1 = double(im1);    im2 = double(im2);
end

if ~exist('patch_wid','var')
    patch_wid = 8; %patch: 8x8
end
Knn_count = size(uv0,3)/2;

% spatial_w = .002; %spatial weight
% spatial_w = spatial_w * ( 2*patch_rad+1 );
% display_flag = 1;

iter_max = 10;

lambda = 3; % coefficients for spatial constraints

[m n c] = size(im1);
m = m - patch_wid +1;
n = n - patch_wid +1;

% if ~exist('method','var')
%     method = 'generalized_charbonnier';
% end

switch lower(func)
    case 'geman_mcclure'
        this.type = @geman_mcclure;
    case 'huber'
        this.type = @huber;
    case 'lorentzian'
        this.type = @lorentzian;
    case 'quadratic'
        this.type = @quadratic;
    case 'tukey'
        this.type = @tukey;
    case 'spline'
        this.type = @spline;
    case 'mixture'
        this.type = @mixture;
    case 'gaussian'
        this.type = @gaussian;
    case 'tdist'
        this.type = @tdist;
    case 'tdist_unnorm'
        this.type = @tdist_unnorm;
    case 'charbonnier'
        this.type = @charbonnier;
    case 'generalized_charbonnier'
        this.type = @generalized_charbonnier;
    otherwise
        error('Invalid robust function type.');
end
this.param = 1;
if strcmp(func, 'generalized_charbonnier')
  this.param = [1e-3, 1];
end;

this.type='gaussian';
this.param=1;

%% compute error map
error_map = zeros(m,n,Knn_count);
for knn_id = 1 : Knn_count
    tmp_uv = uv0( :, :, 2*(knn_id-1)+1 : 2*knn_id );
    for i=1:m
        for j=1:n
            % current location
            tmp_y = i+tmp_uv(i,j,2);
            tmp_x = j+tmp_uv(i,j,1);
            if tmp_x>=1 && tmp_x <= m && tmp_y>=1 && tmp_y<=n%new location is valid
                % compute flow
                tmp_block1 = im1(i:i+patch_wid-1,j:j+patch_wid-1,:);
                tmp_block2 = im2(tmp_y:tmp_y+patch_wid-1,tmp_x:tmp_x+patch_wid-1,:);
                % evaluate error
                tmp_error = evaluate(this, tmp_block1-tmp_block2);
                error_map(i,j,knn_id)=sum(tmp_error(:));
            else
                error_map(i,j,knn_id)=1e3;
            end
        end
    end
end

%% new_uv
%initialize with the first candidate
new_uv = uv0(:,:,1:2); 
uv = new_uv;
tmp_spatial_error = zeros(Knn_count,1);
for iter = 1:iter_max % belief propagation 
    disp('the %d th round of Bp.',iter)
    for i = 1:m
        for j = 1:n
            % data term
            tmp_error_d = reshape(error_map(i,j,:),[Knn_count 1]);
            
            % evaluate spatial constraint error
            tmp_spatial_error(:) = 0;
            if i~=1
                tmp_u = uv0(i,j,1:2:end); tmp_v = uv0(i,j,2:2:end);
                tmp_u1 = uv0(i-1,j,1:2:end); tmp_v1 = uv0(i-1,j,2:2:end);
                tmp_error_s = evaluate(this,tmp_u-tmp_u1) + evaluate(this,tmp_v-tmp_v1);
                tmp_error_s = reshape(tmp_error_s,[Knn_count 1]);
                tmp_spatial_error(:) = tmp_spatial_error(:) + reshape( tmp_error_s, [Knn_count 1]);
            end
            if i~=m
                tmp_u = uv0(i,j,1:2:end); tmp_v = uv0(i,j,2:2:end);
                tmp_u1 = uv0(i+1,j,1:2:end); tmp_v1 = uv0(i+1,j,2:2:end);
                tmp_error_s = evaluate(this,tmp_u-tmp_u1) + evaluate(this,tmp_v-tmp_v1);
                tmp_error_s = reshape(tmp_error_s,[Knn_count 1]);
                tmp_spatial_error(:) = tmp_spatial_error(:) + reshape( tmp_error_s, [Knn_count 1]);
            end
            if j~=1
                tmp_u = uv0(i,j,1:2:end); tmp_v = uv0(i,j,2:2:end);
                tmp_u1 = uv0(i,j-1,1:2:end); tmp_v1 = uv0(i,j-1,2:2:end);
                tmp_error_s = evaluate(this,tmp_u-tmp_u1) + evaluate(this,tmp_v-tmp_v1);
                tmp_error_s = reshape(tmp_error_s,[Knn_count 1]);
                tmp_spatial_error(:) = tmp_spatial_error(:) + reshape( tmp_error_s, [Knn_count 1]);
            end
            if j~=n
                tmp_u = uv0(i,j,1:2:end); tmp_v = uv0(i,j,2:2:end);
                tmp_u1 = uv0(i,j+1,1:2:end); tmp_v1 = uv0(i,j+1,2:2:end);
                tmp_error_s = evaluate(this,tmp_u-tmp_u1) + evaluate(this,tmp_v-tmp_v1);
                tmp_error_s = reshape(tmp_error_s,[Knn_count 1]);
                tmp_spatial_error(:) = tmp_spatial_error(:) + reshape( tmp_error_s, [Knn_count 1]);
            end
            tmp_error = tmp_error_d + lambda * tmp_error_s;
            
            % take the minimum
            [tmp,ind] = min(tmp_error);
            
            % update new_uv
            new_uv(i,j,:) = uv0(i,j,2*(ind-1)+1:2*ind);
        end
    end
    % update uv
    uv = new_uv;
end

end
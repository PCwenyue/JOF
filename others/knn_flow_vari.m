function uv = knn_flow_vari(im1,im2,dominant_uv,top_homo)
% im1, im2: two images
% dominant_uv: 2 x Knn_count potential velocities
% recursively call knn_flow_base to solve dense registration problem

addpath('./flow_code_v2');
addpath(genpath('./flow_code_v2/utils'));
addpath('./gco-v3.0/matlab');
addpath('./gco-v3.0/matlab/bin');
addpath ./qpbo;

%% parameter setting
this = load_of_method('classic+nl');
this.display = 0;
texture = 0;
%this.alp = .5;this.auto_level = 1;
amp_factor = 100;
Knn_count1 = size(dominant_uv,2);
if ~isempty(top_homo)
    Knn_count2 = length(top_homo);
else
    Knn_count2 = 0;
end
%Knn_count = Knn_count1 + Knn_count2;

%% texture-structure decomposition
color_im1 = im1;
if texture      
% Perform ROF structure texture decomposition 
      im1  = structure_texture_decomposition_rof( double(rgb2gray(im1)), 1/8, 100, this.alp);
      im2  = structure_texture_decomposition_rof( double(rgb2gray(im2)), 1/8, 100, this.alp);
      disp('Done with texture-structure decomposition');
else
    im1 = double(rgb2gray(im1));
    im2 = double(rgb2gray(im2));
end;
[m n c] = size(im1);

%% Dominant Velocity Field
if Knn_count1 > 0
    disp('Dominant uv part.');
    %% data term:
    %[x,y]   = meshgrid(1:n,1:m);
    error_map = zeros(Knn_count1,m*n);%error_map = zeros(m,n,Knn_count);
    for i = 1 : Knn_count1
        fprintf('.');
        tmp_uv = dominant_uv(:,i);
        
        [ttmp_uv tmp_error_map]=vary_diff(im1,im2,tmp_uv(1),tmp_uv(2));
        vel(i).u = ttmp_uv(:,:,1);
        vel(i).v = ttmp_uv(:,:,2);
        
        error_map(i,:) = tmp_error_map(:)';
    end
    clear tmp_error_map;
    error_map = round(error_map*amp_factor);
    
    %% spatial term: Total Variation
    cost_matrix = zeros(Knn_count1,Knn_count1);
    for i=1:Knn_count1
        tmp_u = dominant_uv(1:2,i);
        for j=1:Knn_count1
            tmp_v = dominant_uv(1:2,j);
            cost_matrix(i,j) = sum((tmp_u-tmp_v).^2);
        end
    end
    cost_matrix = cost_matrix.^0.4;
    cost_matrix = cost_matrix * this.lambda;
    cost_matrix = round(amp_factor * cost_matrix /10);
    %% Labeling by graph-cut
    
    h = GCO_Create(m*n,Knn_count1); % create handles with mxn pixels to label and Knn_count potential solutions for each pixel
    GCO_SetDataCost(h,error_map); % data term

    % neighbor term
    GCO_SetSmoothCost(h,cost_matrix);   % 

    iy = zeros((m-1)*n,1);
    for i=1:n
        tmp = (i-1)*m+1 : i*m-1;
        iy((i-1)*(m-1)+1:i*(m-1)) = tmp';
    end
    jy = iy + 1;
    ix = 1:m*(n-1);jx = ix + m;
    ifinal = [iy;ix'];jfinal = [jy;jx'];
    [weights_x weights_y] = weighted_derivative_ops_color(m, n, color_im1, 3, .01, 0);
    weights_x = round(weights_x*10); weights_y = round(weights_y*10);
    S = sparse(ifinal,jfinal,[weights_y;weights_x],m*n,m*n);
    GCO_SetNeighbors(h,S);

    % optimization
    tic;    GCO_Expansion(h);    toc;
    uv_label=GCO_GetLabeling(h);
    uv_label=reshape(uv_label,[m n]);

    u = zeros(m,n);
    v = zeros(m,n);
    for i=1:Knn_count1
        tmp_l = find(uv_label==i);
        u(tmp_l) = vel(i).u(tmp_l);
        v(tmp_l) = vel(i).v(tmp_l);
    end
    uv1 = cat(3,u,v);
    GCO_Delete(h);
    clear vel error_map
end

%% Homography part
if Knn_count2 >= 2
    [x,y]   = meshgrid(1:n,1:m);
    error_map = zeros(Knn_count2,m*n);
    for i = 1:Knn_count2
        tmp_homo = top_homo(i).matrix;
        x2 = tmp_homo(1,1)*x + tmp_homo(1,2)*y + tmp_homo(1,3);
        y2 = tmp_homo(2,1)*x + tmp_homo(2,2)*y + tmp_homo(2,3);
        z2 = tmp_homo(3,1)*x + tmp_homo(3,2)*y + tmp_homo(3,3);
        x2 = x2./z2; y2 = y2./z2;

        velocity(i).u = x2-x;    velocity(i).v = y2-y;

        warpIm = zeros(m,n,c);
        for j = 1 : c
            warpIm(:,:,j) = interp2(x,y,im2(:,:,j),x2,y2,'bicubic');
        end

        tmp_error_map = abs(im1-warpIm).^0.8;
        tmp_error_map = inpaintnan(tmp_error_map);
        error_map(i,:) = tmp_error_map(:)';
        fprintf('.');
    end
    clear tmp_error_map;
    error_map = round(error_map*amp_factor);

    %% spatial term: Total Variation
    cost_matrix = zeros(Knn_count2,Knn_count2);
    max_diff = 10;
    cost_matrix(:,:) = max_diff;
    for i=1:Knn_count2
        cost_matrix(i,i) = 0;
    end
    cost_matrix = cost_matrix * this.lambda;
    cost_matrix = round(amp_factor * cost_matrix /10);
    %% graph-cut
    h = GCO_Create(m*n,Knn_count2); % create handles with mxn pixels to label and Knn_count potential solutions for each pixel

    GCO_SetDataCost(h,error_map); % data term

    % neighbor term
    GCO_SetSmoothCost(h,cost_matrix);   % 
    iy = zeros((m-1)*n,1);
    for i=1:n
        tmp = (i-1)*m+1 : i*m-1;
        iy((i-1)*(m-1)+1:i*(m-1)) = tmp';
    end
    jy = iy + 1;
    ix = 1:m*(n-1);jx = ix + m;
    ifinal = [iy;ix'];jfinal = [jy;jx'];
    [weights_x weights_y] = weighted_derivative_ops_color(m, n, color_im1, 3, .01, 0);
    weights_x = round(weights_x*10); weights_y = round(weights_y*10);
    
    S = sparse(ifinal,jfinal,[weights_y;weights_x],m*n,m*n);
    GCO_SetNeighbors(h,S);

    % optimization
    tic;    GCO_Expansion(h);    toc;
    uv_label=GCO_GetLabeling(h);
    uv_label=reshape(uv_label,[m n]);

    u = zeros(m,n);
    v = zeros(m,n);
    for i=1:Knn_count2
        tmp_l = find(uv_label==i);        u(tmp_l) = velocity(i).u(tmp_l);        v(tmp_l) = velocity(i).v(tmp_l);
    end
    uv2 = cat(3,u,v);
    GCO_Delete(h);
    clear velocity
end

%% final fusion by QPBO
if Knn_count1 > 0 && Knn_count2 > 1 % both dominant velocity and homography map exist
    disp('QPBO.');
    uv = qpbo_fusion(im1,im2,uv1,uv2,this);
elseif Knn_count1 ==0 && Knn_count2 > 1
    uv = uv2;
elseif Knn_count1 > 0 && Knn_count2 <= 1
    uv = uv1;
end

end


function [uv error_map]=vary_diff(im1,im2,tmp_u,tmp_v)
[m,n,c] = size(im1);
[x,y]   = meshgrid(1:n,1:m);
uv = [];
count = 1;
for u=-.2:.2:.2
    for v=-.2:.2:.2
        uv = [uv;tmp_u+u tmp_v+v];
    end
end
error_map = zeros(m,n,size(uv,1));

for i=1:9
    tmp_uv = uv(i,:);
    x2 = x + tmp_uv(1); y2 = y+tmp_uv(2);
    warpIm = interp2(x,y,im2(:,:),x2,y2,'bicubic');
    tmp_error_map = abs(im1-warpIm).^0.8;
    if tmp_uv(1)>0,        tmp_error_map(:,end-ceil(tmp_uv(1))+1:end)=repmat(tmp_error_map(:,end-ceil(tmp_uv(1))),[1 ceil(tmp_uv(1))]);
    elseif tmp_uv(1)<0,        tmp_error_map(:,1:-floor(tmp_uv(1)))=repmat(tmp_error_map(:,-floor(tmp_uv(1))+1),[1 -floor(tmp_uv(1))]);
    end
    if tmp_uv(2)>0,        tmp_error_map(end-ceil(tmp_uv(2))+1:end,:)=repmat(tmp_error_map(end-ceil(tmp_uv(2)),:),[ceil(tmp_uv(2)) 1]);
    elseif tmp_uv(2)<0,        tmp_error_map(1:-floor(tmp_uv(2)),:)=repmat(tmp_error_map(-floor(tmp_uv(2))+1,:),[-floor(tmp_uv(2)) 1]);
    end
    error_map(:,:,count) = tmp_error_map;
    count = count +1;
end

[error_map,I] = min(error_map,[],3);

u = zeros(m,n);v = zeros(m,n);
for i=1:size(uv,1)
    tmp_l = find(I==i);        u(tmp_l) = uv(i,1);        v(tmp_l) = uv(i,2);
end

uv = cat(3,u,v);

end
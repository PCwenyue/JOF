function [uv,uv1,uv2] = knn_flow(im1,im2,dominant_uv,top_homo)
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
%texture = 1;
%this.alp = .5;this.auto_level = 1;
this.lambda = this.lambda;
amp_factor = 35;
Knn_count1 = size(dominant_uv,2);
if ~isempty(top_homo)
    Knn_count2 = length(top_homo);
else
    Knn_count2 = 0;
end
%Knn_count = Knn_count1 + Knn_count2;

%% texture-structure decomposition
color_im1 = im1;
[m n c] = size(color_im1);
%if texture      
% Perform ROF structure texture decomposition
% gim1 = double(rgb2gray(im1));gim2 = double(rgb2gray(im2));
% tim1  = structure_texture_decomposition_rof( double(im1), 1/8, 100, this.alp);
% tim2  = structure_texture_decomposition_rof( double(im2), 1/8, 100, this.alp);
disp('Done with texture-structure decomposition');
if c>1
    %im1 = double(rgb2gray(im1));im2 = double(rgb2gray(im2));
    im1 = double(im1);im2 = double(im2);
else
    im1 = double(im1);im2 = double(im2);
end
[m n c] = size(im1);
%% Dominant Velocity Field
if Knn_count1 > 0
    disp('Dominant uv part.');
    %% data term:
    [x,y]   = meshgrid(1:n,1:m);
    error_map = zeros(Knn_count1,m*n);%error_map = zeros(m,n,Knn_count);
    for i = 1 : Knn_count1
        tmp_uv = dominant_uv(:,i);
        
        x2 = x + tmp_uv(1);    y2 = y + tmp_uv(2);

        warpIm = zeros(m,n,c);
        for j = 1 : c
            warpIm(:,:,j) = interp2(x,y,im2(:,:,j),x2,y2,'bicubic');
        end

        tmp_error_map = abs(im1-warpIm).^0.8;    tmp_error_map = mean(tmp_error_map,3);
        %tmp_error_flag = find(isnan(tmp_error_map));    tmp_error_map(tmp_error_flag) = 0;% tmp_error_map .* double(tmp_error_flag);

        if tmp_uv(1)>0,        tmp_error_map(:,end-tmp_uv(1)+1:end)=repmat(tmp_error_map(:,end-tmp_uv(1)),[1 tmp_uv(1)]);
        elseif tmp_uv(1)<0,        tmp_error_map(:,1:-tmp_uv(1))=repmat(tmp_error_map(:,-tmp_uv(1)+1),[1 -tmp_uv(1)]);
        end
        if tmp_uv(2)>0,        tmp_error_map(end-tmp_uv(2)+1:end,:)=repmat(tmp_error_map(end-tmp_uv(2),:),[tmp_uv(2) 1]);
        elseif tmp_uv(2)<0,        tmp_error_map(1:-tmp_uv(2),:)=repmat(tmp_error_map(-tmp_uv(2)+1,:),[-tmp_uv(2) 1]);
        end

        error_map(i,:) = tmp_error_map(:)';
    end
    clear tmp_error_map;
    error_map = round(error_map*100);
    
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
    cost_matrix = cost_matrix * this.lambda; %* size(im1,3);
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
    if c == 3
        [weights_x weights_y] = weighted_derivative_ops_color(m, n, color_im1, 3, .01, 0);
    else
        [weights_x weights_y] = weighted_derivative_ops_grayscale(m, n, color_im1, 3, .01, 0);
    end
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
        u(tmp_l) = dominant_uv(1,i);
        v(tmp_l) = dominant_uv(2,i);
    end
    uv1 = cat(3,u,v);
    GCO_Delete(h);
end

%% Homography part
if Knn_count2 >= 2
    [x,y]   = meshgrid(1:n,1:m);
    error_map = zeros(Knn_count2,m*n);
    for i = 1:Knn_count2
        tmp_homo = top_homo(i).matrix;
        % velocity map generation
        x2 = tmp_homo(1,1)*x + tmp_homo(1,2)*y + tmp_homo(1,3);
        y2 = tmp_homo(2,1)*x + tmp_homo(2,2)*y + tmp_homo(2,3);
        if size(tmp_homo,1) == 3
            z2 = tmp_homo(3,1)*x + tmp_homo(3,2)*y + tmp_homo(3,3);
            x2 = x2./z2; y2 = y2./z2;
        end
        velocity(i).u = x2-x;    velocity(i).v = y2-y;
        
        warpIm = zeros(m,n,c);
        for j = 1 : c
            warpIm(:,:,j) = interp2(x,y,im2(:,:,j),x2,y2,'bicubic');
        end

        tmp_error_map = abs(im1-warpIm).^0.8; tmp_error_map = mean(tmp_error_map,3);
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
    cost_matrix = cost_matrix * this.lambda;% * size(im1,3);
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
    if c == 3
        [weights_x weights_y] = weighted_derivative_ops_color(m, n, color_im1, 3, .01, 0);
    else
        [weights_x weights_y] = weighted_derivative_ops_grayscale(m, n, color_im1, 3, .01, 0);
    end
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
end

%% final fusion by QPBO
if Knn_count1 > 0 && Knn_count2 > 1 % both dominant velocity and homography map exist
    disp('QPBO.');
    uv = qpbo_fusion(im1,im2,color_im1,uv1,uv2,this);
elseif Knn_count1 ==0 && Knn_count2 > 1
    uv = uv2; uv1 = [];
elseif Knn_count1 > 0 && Knn_count2 <= 1
    uv = uv1; uv2 = [];
end

end

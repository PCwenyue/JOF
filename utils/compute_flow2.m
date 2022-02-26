function uv = compute_flow2(this, init, flo_fuse, deqing_flag)
% this: structure contains everything including images
% init: initial flow
% gt: ground truth for evaluation, can be set as empty
if exist('flo_fuse','var') && ~isempty(flo_fuse)   % 1
    fusion_flag = 1; % use fusion
end
% if ~exist('deqing_flag','var')
%     deqing_flag = 0;
% end

  % Frame size
  sz = [size(this.images, 1), size(this.images, 2)];

  % If we have no initialization argument, initialize with all zeros
  if (nargin < 2),    uv = zeros([sz, 2]);
  else    uv = init;
  end
   
  % Preprocess input (gray) sequences
  if this.texture      
      % Perform ROF structure texture decomposition 
      images  = structure_texture_decomposition_rof( this.images, 1/8, 100, this.alp);
      %only the texture part is maintained; the structure part has been discarded  
  elseif this.fc     
      % Laplacian in flowfusion
      f = fspecial('gaussian', [5 5], 1.5);      %f = fspecial('gaussian', [3 3], 1);
      images  = this.images- this.alp*imfilter(this.images, f, 'symmetric');
      
      images  = scale_image(images, 0, 255);
  else
      images  = scale_image(this.images, 0, 255);
  end;
  
%   %%%%%%%%%%%%%%%%%%%%%%%%%改成非对称金字塔%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if 1
%         this.pyramid_levels  =  1 + floor( log(max(sz)/16)...
%             / log(this.pyramid_spacing) );
%         
%         tmp = exp(log(min(sz)/max(sz)*this.pyramid_spacing^(this.pyramid_levels-1))...
%             /(this.pyramid_levels-1) );
%         if sz(1) > sz(2)
%             spacing = [this.pyramid_spacing tmp];
%         else
%             spacing = [tmp this.pyramid_spacing];
%         end
%         
%     end  
%     factor            = sqrt(2);
% if 1
%     fprintf('unequal sampling\n');
%     pyramid_images    = compute_image_pyramid_unequal(images, ...   % 数值金字塔
%         this.pyramid_levels, spacing, factor);
%     
%     % For segmentation purpose
%     org_pyramid_images = compute_image_pyramid_unequal(this.images,... % 彩色金字塔
%         this.pyramid_levels, spacing, factor);
%     org_color_pyramid_images = compute_image_pyramid_unequal(...
%         this.color_images, this.pyramid_levels, spacing, factor);
% end
% For gnc stage 2 to gnc_iters
% smooth_sigma      = sqrt(this.gnc_pyramid_spacing)/factor;
% f                 = fspecial('gaussian',2*round(1.5*smooth_sigma) +1, smooth_sigma);
% gnc_pyramid_images= compute_image_pyramid(images,f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing);

% % For segmentation/weighted median filtering
% org_gnc_pyramid_images = compute_image_pyramid(this.images,...
%     f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing);
% org_color_gnc_pyramid_images = compute_image_pyramid(this.color_images,...
%     f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
  if this.auto_level
      % Automatic determine pyramid level
      this.pyramid_levels  =  1 + floor( log(min(size(images, 1), size(images,2))/16) / log(this.pyramid_spacing) );
  end;
  
  % Construct image pyramid, using filter setting in Bruhn et al in "Lucas/Kanade.." (IJCV2005') page 218
  
  % For gnc stage 1    
  factor            = sqrt(2);  % sqrt(3) worse
  smooth_sigma      = sqrt(this.pyramid_spacing)/factor;   % or sqrt(3) recommended by Manuel Werlberger 
  f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);        
  pyramid_images    = compute_image_pyramid(images, f, this.pyramid_levels, 1/this.pyramid_spacing);


  % For segmentation purpose
  org_pyramid_images = compute_image_pyramid(this.images, f, this.pyramid_levels, 1/this.pyramid_spacing);
  org_color_pyramid_images = compute_image_pyramid(this.color_images, f, this.pyramid_levels, 1/this.pyramid_spacing);

  % For gnc stage 2 to gnc_iters  
  smooth_sigma      = sqrt(this.gnc_pyramid_spacing)/factor;
  f                 = fspecial('gaussian', 2*round(1.5*smooth_sigma) +1, smooth_sigma);        
  gnc_pyramid_images= compute_image_pyramid(images, f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing); 

    
  % For segmentation purpose
  org_gnc_pyramid_images = compute_image_pyramid(this.images, f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing);   
  org_color_gnc_pyramid_images = compute_image_pyramid(this.color_images, f, this.gnc_pyramid_levels, 1/this.gnc_pyramid_spacing);   
  
  %% Levin's modification 2: SIFT based sparse flow propagation
%   if fusion_flag == 1
%       max_flow = max(max(abs(uv(:))),10);
%       if strcmpi(fusion_method,'sift')
%           mymatch = siftmatch(this.images(:,:,1), this.images(:,:,2));
%           delete('tmp.key');    delete('tmp.pgm');
%       else
%           cand_uv = [];
%           for k=1:2 %width 3x3 and 5x5
%               tmp_cand_uv = knn_flow(this.images(:,:,1),this.images(:,:,2),[],max_flow,k);
%               cand_uv =cat(3,cand_uv,tmp_cand_uv);
%           end
%           [m,n,c] = size(cand_uv); knn_count = c/2;
%       end
%    end
  
  for ignc = 1:this.gnc_iters     
  
      if this.display
          disp(['GNC stage: ', num2str(ignc)])
      end
      
      % Update GNC parameters (linearly)
      if this.gnc_iters > 1
          new_alpha  = 1 - (ignc-1) / (this.gnc_iters-1);
          this.alpha = min(this.alpha, new_alpha);
          this.alpha = max(0, this.alpha);          
      end;
      
      if ignc == 1
          pyramid_levels = this.pyramid_levels;          pyramid_spacing = this.pyramid_spacing;
      else
          pyramid_levels = this.gnc_pyramid_levels;          pyramid_spacing = this.gnc_pyramid_spacing;
      end;
      if ignc == 1
          thr = linspace(1e-2,2e-3,pyramid_levels);
      else
          thr = [5e-3 5e-3];
      end
     
      % Iterate through all pyramid levels starting at the top
      for l = pyramid_levels:-1:1
          %% Levin's modification 1: reweighted spatial term, exp(-|\nabla I0|)
          if ignc == 1
              [m n c] = size( org_color_pyramid_images{l} );
              if c == 3,    [weightx weighty] = weighted_derivative_ops_color(m, n, org_color_pyramid_images{l}, 3, 0.01);
              else  [weightx weighty] = weighted_derivative_ops_grayscale(m, n, org_color_pyramid_images{l}, 3, 0.01);
              end
          else
              [m n c] = size( org_color_gnc_pyramid_images{l} );
              if c == 3,    [weightx weighty] = weighted_derivative_ops_color(m, n, org_color_gnc_pyramid_images{l}, 3, 0.01);
              else  [weightx weighty] = weighted_derivative_ops_grayscale(m, n, org_color_gnc_pyramid_images{l}, 3, 0.01);
              end
          end
          
          if this.display
              disp(['-Pyramid level: ', num2str(l)])
          end

          % Generate copy of algorithm with single pyramid level and the
          % appropriate subsampling
          small = this;
          
          if ignc == 1
              nsz               = [size(pyramid_images{l}, 1) size(pyramid_images{l}, 2)];
              small.images      = pyramid_images{l};                            
              small.max_linear  = 1;             % number of linearization performed per warping, 1 OK for quadratic formulation
              im1   = org_pyramid_images{l}(:,:,1);
              small.color_images      = org_color_pyramid_images{l};              
          else
              small.images         = gnc_pyramid_images{l};
              nsz   = [size(gnc_pyramid_images{l}, 1) size(gnc_pyramid_images{l}, 2)];
              im1   = org_gnc_pyramid_images{l}(:,:,1);
              
              small.color_images      = org_color_gnc_pyramid_images{l};
          end;
          
          % Rescale the flow field
          uv        = resample_flow(uv, nsz);  %原始
%            uv = resample_flow_unequal(uv, nsz);
          if exist('gt','var') && ~isempty(gt) && this.display
              tmp_gt = resample_flow(gt, nsz);  [epe ae]=eva_flow2(uv,tmp_gt);
              fprintf('Before optimization: %f %f\n',epe,ae);
          end;
          
          %% Levin's modification2: fusion based initialization
          if fusion_flag == 1 && ignc == this.gnc_iters
              uv2 = resample_flow(flo_fuse,nsz);
              uv = qpbo_fusion(small.images(:,:,1),small.images(:,:,2),small.color_images,uv,uv2,this);
          end

          small.seg = im1;      
          
         % Adaptively determine half window size for the area term         
          small.affine_hsz      = min(4, max(2, ceil(min(nsz)/75)) );
          
          % Run flow method on subsampled images
%           if deqing_flag == 1,             
         uv = compute_flow_base2(small, uv, weightx, weighty);
%           elseif deqing_flag == 0 && ignc ==1 % our fast algorithm
%               uv = compute_flow_base_split_breg_update(small, uv, weightx, weighty, 'ani', 'approx', 1, 5, 1e-2,0);
%           elseif ignc >= 2 % Deqing's cost function with our fast treatment
%               uv = compute_flow_base_split_breg_update(small, uv, weightx, weighty, 'ani', 'approx', 1, 500, thr(l),1);
%           end
%           if this.display && ~isempty(gt)
%               tmp_gt = resample_flow(gt, nsz);  [epe ae]=eva_flow2(uv,tmp_gt);
%               fprintf('After optimzization: %f %f\n',epe,ae);
%               testimg = flowToColor(uv);       imshow(testimg)
%           end
      end
      
%       if this.display
%         fprintf('GNC stage %d finished, %3.2f minutes passed \t \t', ignc, toc/60);
%       end

  end; 
    
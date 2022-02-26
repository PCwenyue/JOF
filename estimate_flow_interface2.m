function uvo = estimate_flow_interface2(im1, im2, method, flo_init, flo_fuse ,params) 
% im1, im2: input images, uint8 rgb/gray, range from 0-255
% method: 'classic+nl'
% flo_init: initialization, if empty set as zeros
% flo_gt: ground truth, can be set as empty
% params: parameters

% Read in arguments
if nargin < 3
    method = 'classic+nl-fast';
end;
if ~exist('flo_init','var')||isempty(flo_init)
     flo_init=zeros([size(im1,1) size(im1,2) 2]);
end;
if ~exist('flo_fuse','var')
    flo_fuse = [];
end
if (~isdeployed)
    addpath('./qpbo');
    addpath('./utils');
    addpath('./flow_code_v2');
    addpath(genpath('./flow_code_v2/utils'));
end

% Load default parameters
ope = load_of_method(method);

if exist('params','var')
    ope = parse_input_parameter(ope, params);    
end;

% Uncomment this line if Error using ==> \  Out of memory. Type HELP MEMORY for your option.
%ope.solver    = 'pcg';  

if size(im1, 3) > 1
    tmp1 = double(rgb2gray(uint8(im1)));
    tmp2 = double(rgb2gray(uint8(im2)));
    ope.images  = cat(length(size(tmp1))+1, tmp1, tmp2);
else    
    if isinteger(im1);
        im1 = double(im1);
        im2 = double(im2);
    end;
    ope.images  = cat(length(size(im1))+1, im1, im2);
end;

% Use color for weighted non-local term
if ~isempty(ope.color_images)    
    if size(im1, 3) > 1        
        % Convert to Lab space       
        im1 = RGB2Lab(im1);          
        for j = 1:size(im1, 3);
            im1(:,:,j) = scale_image(im1(:,:,j), 0, 255);
        end;        
    end;    
    ope.color_images   = im1;
end;

% Compute flow field
%ope.display = 1;
uv  = compute_flow2(ope, flo_init, flo_fuse);

if nargout == 1
    uvo = uv;
end;

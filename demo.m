
addpath(genpath('./csh')); 
addpath ./flow-code-matlab;
addpath ./utils;
I0 = imread('training\clean2\alley_2\frame_0032.png');
I1 = imread('training\clean2\alley_2\frame_0033.png');
%% direct CSH
csh_wid = 8;
[uv1,flag1] = CSH_nn_flow_base(I0,I1,[],csh_wid,1);
[histmat_csh top_uv2] = Dominant_Offset(uv1,flag1);

%% top homography computation by SIFT
top_homo = Dominant_Transform(uv1);
    
%% motion segmentation
[uv,uv1,uv2] = knn_flow(I0,I1,top_uv2,top_homo);

%% continuous refinement
uvo = estimate_flow_interface2(I0,I1, 'classic+nl', [], uv);
toc;
flow = readFlowFile('training\flow\alley_2\frame_0032.flo');
tu = flow(:,:,1);
tv = flow(:,:,2);
                     
UNKNOWN_FLOW_THRESH = 1e9; 
tu (tu>UNKNOWN_FLOW_THRESH) = NaN;
tv (tv>UNKNOWN_FLOW_THRESH) = NaN;
imshow(uint8(flowToColor(uvo))); 


if sum(~isnan(tu(:))) > 1

    [aae stdae aepe] = flowAngErr(tu, tv, uvo(:,:,1), uvo(:,:,2), 0); % ignore 0 boundary pixels
    fprintf('\nAAE %3.3f average EPE %3.3f \n', aae, aepe);
    
end;
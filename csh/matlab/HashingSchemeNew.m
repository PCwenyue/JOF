function [AnnA2B,AnnB2A,...     % results, each is 3 times the size of A or B
    innerTimes,...              % profiling
    annErrorsA,annErrorsB]...   % of size A*(iterations+1) - hold errors after each iteration (and before beginning)
    = HashingSchemeNew(rgbA,rgbB,k,calcBnn,parameters,visualizationType,mask,patch_mode)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % parameters.width                                -> Kernel width (for non-descriptor mode)
% % parameters.descriptor_params.descriptor_mode    -> Boolean saying is this patch or descriptor mode
% % parameters.descriptor_params.hA -> Height of image A descriptors
% % parameters.descriptor_params.wA -> Width of image A descriptors
% % parameters.descriptor_params.hB -> Height of image B descriptors
% % parameters.descriptor_params.wB -> Width of image B descriptors

% % parameters.dist
% % parameters.maxKernels -> Maximum number of kernels for calculations
%%  parameters.useTicsInside -> Do profiling
% % parameters.TLboundary -> Top left boundary
% % parameters.BRboundary -> Bottom right boundary
% % parameters.br_boundary_to_ignore -> Bottom right boudaries to ignore in the mapping calculation
% % parameters.numHashTables -> Number of hash tables to use
% % parameters.numHashs -> Width of hash table (should be 2)
% % parameters.insideInfo -> Print debug info


%% A] extracting input paramters

if (nargin < 5)
    error('HashingSchemeNew error: At least 5 input paraemters are needed');
end
width = parameters.width; %8
descriptor_mode = parameters.descriptor_params.descriptor_mode; %0
maxKernels = parameters.maxKernels; %27
useTics = parameters.useTicsInside; %0
TLboundary = parameters.TLboundary; %16
BRboundary = parameters.BRboundary; %8
br_boundary_to_ignore = parameters.br_boundary_to_ignore; %8
numTables = parameters.numHashTables; %5
numHashs = parameters.numHashs; %2
insideInfo = parameters.insideInfo;
CalcErrorImages = parameters.DebugCalcErrors;
dist = parameters.dist; % this is for CSH that keeps 'dist' pixels in x or y away from identity


% timers
if (useTics)
    innerTimes.PreProcessTime = 0;
    innerTimes.GCKTime = 0;
    innerTimes.ModifyBoundariesTime = 0;
    innerTimes.MemoryConvertTime = 0;
    innerTimes.CodesTime = 0;
    innerTimes.HashingTime = 0;
    innerTimes.CandidateCreation = 0;
    innerTimes.IterationTimes = zeros(1,numTables);
    innerTimes.PostProcessTime = 0;
end

if (useTics)
    PreProcessTime = tic;
end

rotation_invariant = 0;

if (~descriptor_mode)
    if (~patch_mode)
        % sizes
        [hA,wA,dA] = size(rgbA);
        [hB,wB,dB] = size(rgbB);
        if (ndims(rgbA)== 3)
            A = rgb2ycbcr(rgbA);
        else
            A = rgbA; % Gray level image
        end
        
        if (ndims(rgbB) == 3)
            B = rgb2ycbcr(rgbB);
        else
            B = rgbB; % Gray level image
        end
    else
        A = rgbA; B = rgbB; % the input A has size: size(A) = [width^2,(hA-[width-1])*(wA-[width-1])*dA]
        hA = parameters.descriptor_params.hA;
        wA = parameters.descriptor_params.wA;
        hB = parameters.descriptor_params.hB;
        wB = parameters.descriptor_params.wB;
        dA = parameters.descriptor_params.dA;
        dB = dA;
        rotation_invariant = parameters.descriptor_params.rotation_invariant;
        if (rotation_invariant && isequal(A,B))
            dist = 1;
        end
    end
else % descriptor mode
    if ((ndims(rgbA) ~= 2) || (ndims(rgbB) ~= 2))
        error('Input data must be on the type of (Descripotr_Size*(W*H))');
    end
    
    if (~strcmp(class(rgbA),'single'))
        error('Input data A must be of type single');
    end
    if (~strcmp(class(rgbB),'single'))
        error('Input data B must be of type single');
    end
    
    A = rgbA; B = rgbB;
    hA = parameters.descriptor_params.hA;
    wA = parameters.descriptor_params.wA;
    hB = parameters.descriptor_params.hB;
    wB = parameters.descriptor_params.wB;
    dA = 1; dB = 1;
    width = 1;
end

if (dA ~= dB)
    error ('Image color channels must be the same: both RGB or both gray level images');
end

ColorChannels = dA;

% long term initializations
sizA = [hA,wA]; sizB = [hB,wB];


visualizeCandidates = 0;
visualizeHashTables = 0;

% Set wheather to create Walsh Hadamard projections in memory at the
% optimal way for algorithm run time(NumKernels*Img_Size matrix) or in the
% logical view (Number of images which are the projection results, as the number of WH kernels)

visualizeFilterDistributions = 0;
visualizeDenoising = false;

if (exist('visualizationType','var'))
    if (strcmp(visualizationType,'hashTableView'))
        visualizeHashTables = true;
        figure;
    end
    if (strcmp(visualizationType,'filterDistributions'))
        visualizeFilterDistributions = true;
    end
    if (strcmp(visualizationType,'candidatesView'))
        visualizeCandidates = 1;
        figure;
    end
    if (strcmp(visualizationType,'denoising'))
        visualizeDenoising = 1;
    end
end

if (~exist('mask','var'))
    mask = [];
else
    if (~isempty(mask))
        [hm,wm,dm] = size(mask);
        
            if (hm ~= hB) || (wm ~= wB) | (dm ~= 1)
                error 'mask image must have the same dimensions like target image (iamge B)';
            end
        if (calcBnn)
            if (hm ~= hA) || (wm ~= wA) | (dm ~= 1)
                error 'mask image must have the same dimensions like target image (iamge A)';
            end
        end
        % Convert the original mask which is in "patches" to pixel mask,
        % to do that we should add the pixels in the size of (width -1 * width - 1) in the top left corener of the mask
        % e.g. mask is 0 0 0 0		output mask (width = 2): 0 1 1 0
        %			   0 0 1 0								 0 1 1 0
        %              0 0 0 0								 0 0 0 0
        %			   0 0 0 0								 0 0 0 0
        %
        
        
        if (0) % GUY: Optional only
            se = zeros(width*2-1 , width*2-1);
            se(1:width,1:width) = 1;
            mask = imdilate(mask,se);
        end
    end
end
%% B] PARAMETERS/INITIALIZTIONS 1 -Internal (only for statistics)

% Collect statistics for analysis. Two kinds of information:
% 1) General statistics on number and quality of candidates: a histogram of the number of
%       candidates of each quality (=depth=numFilters)
% 2) Image (per pixel) statistics:
%       a) number of candidates image
%       b) quality of candidates image (average of depths / sum of depths)
%%% data structures for this:
numKernelsPerAlternative = zeros(1,numTables);
bitCountPerTable = zeros(1,numTables);
%%%

% for debug - follow error improvements after each table
if (CalcErrorImages)
    annErrorsA = zeros(1,numTables+1); % since I also save the error before the first iteration
    if (calcBnn)
        annErrorsB = zeros(1,numTables+1); % since I also save the error before the first iteration
    else
        annErrorsB = [];
    end
end

%% C] PARAMETERS/INITIALIZTIONS 2 - Actual things

% depending on the patch width - how many kernels (maximum) do we want to compute

if (descriptor_mode == 0)
    SupportedWidth = [2 4 8 16 32];
    
    if (isempty(find(width == SupportedWidth, 1)))
        SupportedWidthStr = [];
        for L = 1 : length(SupportedWidth)
            SupportedWidthStr = [SupportedWidthStr num2str(SupportedWidth(L)) ' '];
        end
        error(['The input width ' num2str(width) ' is not supported, only width size: ' SupportedWidthStr 'are supported']);
    end
    
    % maxBits = the number of bits in the code (originally was: % 'ceil(log2(hashFact*max(hA*wA,hB*wB)));' )
    
    switch(width)
        case 2
            % GUY: Just for now - to fill after research
            maxKernels = 2*2;
            maxBits = 15;
        case 4
            maxKernels = 3*3;
            maxBits = 17;
        case 8
            maxKernels = 5*5;
            maxBits = 18;
        case 16
            maxKernels = 7*7;
            maxBits = 18;
        case 32
            maxKernels = 9*9;
            maxBits = 20;
    end
    
    maxKernels = maxKernels * ColorChannels;
else % Descriptors mode
    width = 8; % Width for descriptor mode
    [Descriptor_Width_A NumProjections] = size(A);
    if (Descriptor_Width_A  < 50)
        maxKernels = ceil(Descriptor_Width_A/3);
        maxBits = 16;
    else
        if (Descriptor_Width_A  < 100)
            maxKernels = ceil(Descriptor_Width_A/2);
            maxBits = 18;
        else
            maxKernels = ceil(0.7*Descriptor_Width_A);
            maxBits = 24;
        end
    end
end
run_C_code = 1;

% prepare result matrices
d_mapping = 2; % Two dimensions to the image: x and y.
if (insideInfo)
    d_mapping = d_mapping + 1; % Add the error dimnsion in case of debug
end
AnnA2B = ones(hA,wA,d_mapping,k,'int32');
if (calcBnn)
    AnnB2A = ones(hB,wB,d_mapping,k,'int32');
else
    AnnB2A = [];
end
% Choose type
if (descriptor_mode == 0)
    if (width <= 8)
        classType = 'int16';
    else
        classType = 'int32';
    end
else
    classType = 'int16';%'int32';%
end
% explanations of the following variables: (symmetric for A and B)

% - nBestMapping32uA: of size like A, holds the current best found mapping which is by a FLAT index into B, that runs column after column
% - bestErrorsNewA  : of size like A, holds the approximated errors (GCK errors, not SSD errors) of the current best mapping

% alocate memory
nBestMapping32uA = zeros(hA,wA,k,'uint32'); % current best mapping
bestErrorsNewA = zeros(hA,wA,k,'uint32'); % current (GCK)-errors (approximation of SSD error)

if (useTics)
    innerTimes.PreProcessTime = toc(PreProcessTime);
end

% this calculates all filter images (returned in INT16 format)
if (useTics)
    GCK_time =tic;
end
% GRTIME = tic;
% Output PCA projections is in class signal
if (patch_mode) % here [hA,wA] is the size of the (width-1)-BR-padded valid patches image and therefore, each currFiltImgs_A is of size [hA,wA]
    % the input A has size: size(A) = [width^2,(hA-[width-1])*(wA-[width-1])*dA]
    debugDisp(insideInfo, 'Getting Walsh Hadamard GCK projections for patch mode...');
    [currFiltImgs_A,currFiltImgs_B,nSequencyOrder16u,nSequencyLevels16u ,WHK_with_Cb_Cr] = ...
        GetResultsOfKernelApplication(A,B,TLboundary,BRboundary,width,classType,maxKernels,hA,wA,dA,hB,wB);
elseif (descriptor_mode)
    debugDisp(insideInfo, 'getting principle component analisys results...');
    [currFiltImgs_A,currFiltImgs_B,PCA_A,PCA_B,nSequencyOrder16u,nSequencyLevels16u,MaxDescriptorIntVal] = ...
        GetDescriptorPCA(A,B,hA,wA,hB,wB,TLboundary,BRboundary,classType,maxKernels);
    WHK_with_Cb_Cr = []; % to maintain compatibility
else
    debugDisp(insideInfo, 'Getting Walsh Hadamard GCK projections ');
    [currFiltImgs_A,currFiltImgs_B,nSequencyOrder16u,nSequencyLevels16u ,WHK_with_Cb_Cr] = ...
        GetResultsOfKernelApplication(A,B,TLboundary,BRboundary,width,classType,maxKernels);
end
% GR = toc(GRTIME)
%   save 'before_P.mat'
% keyboard
if (k > 1)
    % If we work on the same image, use enrichment mode for candidate
    % ranknig: use the fact that both images are identical to get more
    % candidates
    if (isequal(rgbA , rgbB))
        KNN_enrichment_mode = 1;
    else
        KNN_enrichment_mode = 0;
    end
end

debugDisp(insideInfo, '> done');

if (useTics)
    innerTimes.GCKTime = toc(GCK_time);
end

% initialize nBestMapping32uA with a random mapping (for each k)
GckErrorTime  = 0;
for i = 1 : k
    if (useTics)
        loop_t = tic;
    end
    bestMappingNewA = (i+ max(hA,hB) + width + 1) * ones(sizA,'uint32');
    nBestMapping32uA(:,:,i) = uint32(bestMappingNewA)-1; % C syntax
    % initial errors
    [bestMappingSubbedA(:,:,2),bestMappingSubbedA(:,:,1)] = ind2sub(sizB,bestMappingNewA);
    if (CalcErrorImages)
        [GCKrms1,annErrorsA(1)] = GetErrorMeanAndImage(bestMappingSubbedA,hB,wB,hA,wA,br_boundary_to_ignore,rgbA,rgbB,width,descriptor_mode);
    end
    bestErrorsNewA(:,:,i) = (2^32-1); % Max int, force candidate check
    if (useTics)
        GckErrorTime = GckErrorTime  + toc(loop_t);
    end
end

% repeat the same for B
if (calcBnn)
    if (useTics)
        loop_t = tic;
    end
    nBestMapping32uB = zeros(hB,wB,k,'uint32');
    bestErrorsNewB = zeros(hB,wB,k,'uint32');
    for i = 1 : k
        bestMappingNewB = (i+ max(hB,hA) + width + 1) * ones(sizB,'uint32');
        nBestMapping32uB(:,:,i) = uint32(bestMappingNewB)-1; % C syntax
        % initial errors
        [bestMappingSubbedB(:,:,2),bestMappingSubbedB(:,:,1)] = ind2sub(sizA,bestMappingNewB);
        if (insideInfo)
            [GCKrms1,annErrorsB(1)] = GetErrorMeanAndImage(bestMappingSubbedB,hA,wA,hB,wB,br_boundary_to_ignore,rgbB,rgbA,width,descriptor_mode);
        end
        
        bestErrorsNewB(:,:,i) = (2^32-1); % Max int, force candidate check
        if (useTics)
            GckErrorTime = GckErrorTime  + toc(loop_t);
        end
    end
end

if (useTics)
    innerTimes.GCKTime = innerTimes.GCKTime + GckErrorTime;
    debugShow(insideInfo, 'innerTimes.GCKTime',innerTimes.GCKTime);
end

% This is done to compensate, adding offset value in Build walsh hadamard
% kernels
if (useTics)
    boundaries_time = tic;
end

if (descriptor_mode)
    if (strcmp(classType,'int16'))
        currFiltImgs_A = int16(currFiltImgs_A);
        currFiltImgs_B = int16(currFiltImgs_B);
    else
        currFiltImgs_A = int32(currFiltImgs_A);
        currFiltImgs_B = int32(currFiltImgs_B);
    end
end

%% D] CALCULATING WH KERNELS - all done in advance (to be kept in memory)

% Trick - create big penalty for pixels at image margin (to avoid memory exceptions in propagation)
% Large val should be the maximum possible value of patch projection, which
% is the DC part: SUM(PATH_GRAY_LEVELS). - Since there are 256 gray levels
% (2^8), we get the following formula
if (descriptor_mode == 0)
    largeVal = 256*width*width-1; % This is max positive values range -1(e.g.: in 8 kernel bit, range is 8*8*256=16484, and largeval =  16383
else
    largeVal = MaxDescriptorIntVal;
end

% Prevent overflow in next algorithm stages cause projections are powered
% by 2 for error calculation (error is 32 bit unsigned int).
% So this is how we take care that error would not overflow
MaxVal = 2^16-1;
if (largeVal > MaxVal)
    largeVal = MaxVal;
end

% this first call is to support properly the B2A direction of the bidirectional mode
[currFiltImgs_A_B2A, currFiltImgs_B_B2A] = ModifyBoundaries(currFiltImgs_A, currFiltImgs_B,largeVal,1,hA,hB,width,descriptor_mode,mask);
[currFiltImgs_A, currFiltImgs_B]         = ModifyBoundaries(currFiltImgs_A, currFiltImgs_B,largeVal,0,hA,hB,width,descriptor_mode,mask);

if (useTics)
    innerTimes.ModifyBoundariesTime = innerTimes.ModifyBoundariesTime + toc(boundaries_time);
end

%% E] Build walsh hadamard codes
% description: Build codes by using the WH projections calculated before
% and randmoally selecting a shift

debugDisp(insideInfo, '-----------  Build codes  --------------------------------------------------------------');

useWH_codes = 1; % Initial settings
if (descriptor_mode == 1)
    usePCA_codes = 1;
    useWH_codes = 0; % Do not use Walsh hadamard codes, rather use PCA codes
end
code_A = cell(1,numTables);
code_B = cell(1,numTables);

CodesAddition = zeros(1,numTables);

for tableInd = 1 : numTables % (= iterations)
    
    if (useTics)
        tic
    end
    
    if (useWH_codes)
        % here I make a compromise (this is perfectly suitable only for A2B and not perfect for B2A)
        [code_A_Img,code_B_Img,bitCountPerTable(tableInd),numKernelsPerAlternative(tableInd)] = ...
            BuildWalshHadamardCodes(tableInd,width,currFiltImgs_A, currFiltImgs_B, ...
            visualizeFilterDistributions,numTables,hA,wA,hB,wB,insideInfo,classType , ColorChannels , WHK_with_Cb_Cr);
        
    else
        if (usePCA_codes)
            [code_A_Img,code_B_Img,bitCountPerTable(tableInd),numKernelsPerAlternative(tableInd)] = ...
                Build_PCA_Codes(tableInd,currFiltImgs_A, currFiltImgs_B,visualizeFilterDistributions,numTables,hA,wA,hB,wB,...
                MaxDescriptorIntVal,Descriptor_Width_A,classType);
        else
            [code_A_Img,code_B_Img,bitCountPerTable(tableInd),numKernelsPerAlternative(tableInd)] = ...
                BuildLSHCodes(A,B,tableInd,width,visualizeFilterDistributions,numTables,hA,wA,hB,wB,calcBnn);
        end
    end
    
    code_A{tableInd}=code_A_Img;
    code_B{tableInd}=code_B_Img;
    
    if (useTics)
        CodesAddition(tableInd) = toc;
        innerTimes.CodesTime =  innerTimes.CodesTime + CodesAddition(tableInd);
    end
end

debugDisp(insideInfo, '--');

%% F] Memory convert
% description: Replace the filtering results images in memory with a more
% optimal arrangement in memory. This dramatically decreases candidate
% ranking run time.

debugDisp(insideInfo, '-----------  Memory convert --------------------------------------------------------------');

if (useTics)
    tic
end

if (~descriptor_mode)
    Num_Filters = length(currFiltImgs_A);
    [Filter_A_H Filter_A_W] = size(currFiltImgs_A{1});
    % Allocate memory in the opposite way  to the one needed ([W*H , nKernels]),
    % THIS IS NOT A BUG, but in order to save time. doing so and making
    % transpose at the end is much faster than straight copy.
    CF_Vectors_A = zeros(Filter_A_H * Filter_A_W , Num_Filters, classType);
    for iteration = 1 : Num_Filters
        CF_Vectors_A(:,iteration) = currFiltImgs_A{nSequencyOrder16u(iteration)}(:);
    end
    
    % Transpose the memory block to get the desired memory organization
    CF_Vectors_A = CF_Vectors_A';
    
    if (visualizeCandidates == 0 && ~visualizeDenoising)
        clear currFiltImgs_A; % Clear previous representation of filter projections on image A
    end;
    
    
    [Filter_B_H Filter_B_W] = size(currFiltImgs_B{1});
    CF_Vectors_B = zeros(Filter_B_H * Filter_B_W , Num_Filters , classType);
    for iteration = 1 : Num_Filters
        CF_Vectors_B(:,iteration) = currFiltImgs_B{nSequencyOrder16u(iteration)}(:);
    end
    
    % Transpose the memory block to get the desired memory organization
    CF_Vectors_B = CF_Vectors_B';
    
    if (visualizeCandidates == 0 && ~visualizeDenoising)
        clear currFiltImgs_B; % Clear previous representation of filter projections on image B
    end;
    
    if (calcBnn) % create a second array of vectors (modify-boundaries was different in this case)
        CF_Vectors_A_B2A = zeros(Filter_A_H * Filter_A_W , Num_Filters, classType);
        for iteration = 1 : Num_Filters
            CF_Vectors_A_B2A(:,iteration) = currFiltImgs_A_B2A{nSequencyOrder16u(iteration)}(:);
        end
        
    % Transpose the memory block to get the desired memory organization
        CF_Vectors_A_B2A = CF_Vectors_A_B2A';
        
        if (visualizeCandidates == 0 && ~visualizeDenoising)
            clear currFiltImgs_A_B2A; % Clear previous representation of filter projections on image A
        end;        
        
        CF_Vectors_B_B2A = zeros(Filter_B_H * Filter_B_W , Num_Filters , classType);
        for iteration = 1 : Num_Filters
            CF_Vectors_B_B2A(:,iteration) = currFiltImgs_B_B2A{nSequencyOrder16u(iteration)}(:);
        end
        
        % Transpose the memory block to get the desired memory organization
        CF_Vectors_B_B2A = CF_Vectors_B_B2A';
        
        if (visualizeCandidates == 0 && ~visualizeDenoising)
            clear currFiltImgs_B_B2A; % Clear previous representation of filter projections on image B
        end;
        
    end
    
else %(~descriptor_mode)
    CF_Vectors_A = reshape(currFiltImgs_A,hA*wA,maxKernels)';
    clear currFiltImgs_A;
    CF_Vectors_B = reshape(currFiltImgs_B,hB*wB,maxKernels)';
    clear currFiltImgs_B;
    if (calcBnn)
        CF_Vectors_A_B2A = reshape(currFiltImgs_A_B2A,hA*wA,maxKernels)';
        clear currFiltImgs_A_B2A;
        CF_Vectors_B_B2A = reshape(currFiltImgs_B_B2A,hB*wB,maxKernels)';
        clear currFiltImgs_B_B2A;        
    end
end

if (k > 1)
    if (calcBnn)
        nBestMapping32uB = reshape(nBestMapping32uB, hB*wB,k)';
        bestErrorsNewB = reshape(bestErrorsNewB, hB*wB,k)';
    else
        nBestMapping32uA = reshape(nBestMapping32uA, hA*wA,k)';
        bestErrorsNewA = reshape(bestErrorsNewA, hA*wA,k)';
    end
end

if (useTics)
    innerTimes.MemoryConvertTime = toc;
end

debugDisp(insideInfo, '--');


debugDisp(insideInfo, '-----------  HUGE loop --------------------------------------------------------------');


%% G] HUGE (MAIN) loop
% description: We have at hand all the filters organized in an optimal way
% in the memory for the ranknig stage: matrix of [number_of_filter *  (Width* Hegiht)]
% We construct one hash table at a time, from the previously aved code images.
%Then, we create candidates for each pixel (patch) and  choose the best candidate
% (comparing also the best candidate from the previous round.

MaxHashTableBits = max(bitCountPerTable);

if (MaxHashTableBits < maxBits)
    maxBits = MaxHashTableBits;
end

for tableInd = 1 : numTables % (= iterations)
    % 2] Hash into tables
    %    ----------------
    if (useTics)
        tic
    end
    if (tableInd==1)
        hTables_A = [];hTables_B = [];
    end
    debugDisp(insideInfo, '== Before HashCodesIntoTables');
    [hTables_A,hTables_B,indices_A,indices_B,indices_B_flat,indices_B_offset] = ...
        HashCodesIntoTables(hA,wA,hB,wB,maxBits,bitCountPerTable,tableInd,numHashs,...
        code_A{tableInd},code_B{tableInd},hTables_A,hTables_B,insideInfo);
    debugDisp(insideInfo, '== After HashCodesIntoTables');
    if (useTics)
        HashingAddition = toc;
        innerTimes.HashingTime = innerTimes.HashingTime + HashingAddition;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % visualizations
    if ( visualizeFilterDistributions)
        continue
    end
    
    if (visualizeCandidates == 0 && ~visualizeDenoising)
        currFiltImgs_A = []; currFiltImgs_B = [];
    end
    
    if (insideInfo)
        HashingSchemeNewVisualizations(visualizeHashTables,visualizeCandidates,...
            indices_A,numHashs,rgbA,rgbB,tableInd,code_A{tableInd},hTables_A,code_B{tableInd},hTables_B,...
            nBestMapping32uA(:,:,1),indices_B_flat,indices_B_offset, maxKernels,currFiltImgs_A,currFiltImgs_B);
    end
    %%%%%%%%%%%%%%%%%%%%%%%
    
    % 3] Candidate creation and ranking
    %    ------------------------------
    
    if (useTics)
        tic
    end
    
    if (run_C_code)
        
        nHashsNoe16u = uint16(numHashs);
        
        % this indicates whether we are going upwards or downwards
        bDownwardFlag  = logical(mod(tableInd,2) == 1);
        
        if (k > 1)
            if strcmp(classType,'int16')
                if (KNN_enrichment_mode)
                    alg = @GCKCodebookPropagation16s_KNN_same_image_B;
                else
                    alg = @GCKCodebookPropagation16s_KNN_B;
                end
            else
                if (KNN_enrichment_mode)
                    alg = @GCKCodebookPropagation32s_KNN_same_image_B;
                else
                    alg = @GCKCodebookPropagation32s_KNN_B;
                end
            end
            alg(bestErrorsNewA, nBestMapping32uA, CF_Vectors_A, CF_Vectors_B, indices_A, indices_B, hTables_A, hTables_B, ...
                nHashsNoe16u, nSequencyOrder16u, nSequencyLevels16u, bDownwardFlag, uint32(k));
        else % (k>1)
            if (visualizeDenoising)
                GCKCodebookPropagation16s_1_2_newA_FAR(bestErrorsNewA, nBestMapping32uA, currFiltImgs_A, currFiltImgs_B, indices_A, indices_B, hTables_A, hTables_B, ...
                    nHashsNoe16u, nSequencyOrder16u, nSequencyLevels16u, bDownwardFlag,int16(dist));
            else
                if  (rotation_invariant && dist > 0 && strcmp(classType,'int16'))
                    GCKCodebookPropagation16s_FilterVectors_Rot_Inv_with_dist(...
                        bestErrorsNewA, nBestMapping32uA, CF_Vectors_A, CF_Vectors_B, indices_A, indices_B, hTables_A, hTables_B, ...
                        nHashsNoe16u, nSequencyOrder16u, nSequencyLevels16u, bDownwardFlag,int16(dist));
                else
                    if strcmp(classType,'int16')
                        if (rotation_invariant)
                            alg = @GCKCodebookPropagation16s_FilterVectors_Rot_Inv;
                        else
                            alg = @GCKCodebookPropagation16s_FilterVectors;
                        end
                    else
                        alg = @GCKCodebookPropagation32s_FilterVectors;
                    end
                    % save uniWS.mat
                    alg(bestErrorsNewA, nBestMapping32uA, CF_Vectors_A, CF_Vectors_B, indices_A, indices_B, hTables_A, hTables_B, ...
                        nHashsNoe16u, nSequencyOrder16u, nSequencyLevels16u, bDownwardFlag);
                end
                % save biWS.mat
                if (calcBnn)
                    alg(bestErrorsNewB, nBestMapping32uB, CF_Vectors_B_B2A, CF_Vectors_A_B2A, indices_B, indices_A, hTables_B, hTables_A, ...
                        nHashsNoe16u, nSequencyOrder16u, nSequencyLevels16u, bDownwardFlag);
                end
            end
        end
        if (useTics)
            CandidatesAddition = toc;
            innerTimes.CandidateCreation = innerTimes.CandidateCreation + CandidatesAddition;
            debugShow(insideInfo, '== Times: [Codes, Hashing,Candidates]',[CodesAddition(tableInd) HashingAddition CandidatesAddition]);
            innerTimes.IterationTimes(tableInd) = CodesAddition(tableInd) + HashingAddition + CandidatesAddition;
        end
        
        % debug - collect new errors after this iteration (for k case - take first matrix)
        if (CalcErrorImages)
            
            if (k > 1) % Convert error and mapping to images shape in order to compute GCK error for debugging
                nBestMapping32uA = reshape(nBestMapping32uA', hA,wA,k);
                if (calcBnn)
                    nBestMapping32uB = reshape(nBestMapping32uB', hB,wB,k);
                end
            end
            
            for mapInd = 1 : k
                [bestMappingSubbedA(:,:,2),bestMappingSubbedA(:,:,1)] = ind2sub(sizB,nBestMapping32uA(:,:,mapInd)+1);% C-syntax
                if (descriptor_mode == 0)
                    [GCKrms1,mexError] = GetErrorMeanAndImage(bestMappingSubbedA,hB,wB,hA,wA,br_boundary_to_ignore,rgbA,rgbB,width,0);
                else
                    [GCKrms1,mexError]= GetErrorMeanAndImage(bestMappingSubbedA,hB,wB,hA,wA,br_boundary_to_ignore,A,B,width,1);
                end
                annErrorsA(tableInd+1) = annErrorsA(tableInd+1) + mexError;
            end
            annErrorsA(tableInd+1) = annErrorsA(tableInd+1) ./k;
            debugShow(insideInfo, '============== [annErrorsAvg]',[annErrorsA(tableInd+1)]);
            
            if ((k > 1) && (tableInd < numTables)) % Convert mapping back to format that algo should receive
                if (calcBnn)
                    nBestMapping32uB = reshape(nBestMapping32uB, hB*wB,k)';
                else
                    nBestMapping32uA = reshape(nBestMapping32uA, hA*wA,k)';
                end
            end
            
            debugShow(insideInfo, 'mexError',mexError);
            
            if (calcBnn)
                for mapInd = 1 : k
                    [bestMappingSubbedB(:,:,2),bestMappingSubbedB(:,:,1)] = ind2sub(sizA,nBestMapping32uB(:,:,mapInd)+1);% C-syntax
                    if (descriptor_mode == 0)
                        [GCKrms1,mexError] = GetErrorMeanAndImage(bestMappingSubbedB,hA,wA,hB,wB,br_boundary_to_ignore,rgbB,rgbA,width,0);
                    else
                        [GCKrms1,mexError] = GetErrorMeanAndImage(bestMappingSubbedB,hA,wA,hB,wB,br_boundary_to_ignore,B,A,width,1);
                    end
                    annErrorsB(tableInd+1) = annErrorsB(tableInd+1) + mexError;
                end
                annErrorsB(tableInd+1) = annErrorsB(tableInd+1) ./k;
                debugShow(insideInfo, 'mexError',mexError);
            end
        end
        
        
        debugDisp(insideInfo, ['>>>>>>>>>                           done                       >>>>>>>' char(10)] );
        
        
    else % run c-code
        
        % num Candidates calculus:
        % numHashs (type 1) + 2*numHashs (type 2) + 2*numHashs (type 3) +
        % 2*numHashs (type 4) + 2*numHashs^2 (type 5) + numHashs (type A)
        % We leave out types 3,4, 5 for the moment -> we get: 4*numHashs +  1 for current best
        numCandidates =  4*numHashs+1 + 2;% + 2 for types 2
        
        if (mod(tableInd,2) == 1) % odd table - we run downwards
            for j = 2 : wA
                if (mod(j,100) == 2)
                    debugDisp(insideInfo, ['Forw:  col ' num2str(j) ' out of ' num2str(wA)]);
                end
                
                for i = 2 : hA % here starting from 2 becuase of propagation
                    [candidates,bestInd,bestMappingNewA] = ...
                        DownFunctionMatlab(numCandidates,bestMappingNewA,i,j,...
                        indices_A,hTables_A,hTables_B,...
                        hA,wA,hB,wB,...
                        indices_B_flat,indices_B_offset,...
                        numHashs,maxKernels,currFiltImgs_A,currFiltImgs_B);
                end
            end
        else  %  of  if (mod(tableInd,2) == 1) % odd table
            % so now we run upwards
            for j = wA-1 : -1 : 1
                if (mod(j,100) == 2)
                    debugDisp(insideInfo, ['Back:  col ' num2str(j) ' out of ' num2str(wA)]);
                end
                
                for i = hA-1 : -1 : 1 % here starting from hA-1 becuase of propagation
                    [candidates,bestInd,bestMappingNewA] = ...
                        UpFunctionMatlab(numCandidates,bestMappingNewA,i,j,...
                        indices_A,hTables_A,hTables_B,...
                        hA,wA,hB,wB,...
                        indices_B_flat,indices_B_offset,...
                        numHashs,maxKernels,currFiltImgs_A,currFiltImgs_B);
                    
                end
            end
        end
        
        % debug - collect new errors after this iteration
        [bestMappingSubbedA(:,:,2),bestMappingSubbedA(:,:,1)] = ind2sub(sizB,bestMappingNewA);
        if (descriptor_mode == 0)
            [GCKrms1,regularError] = GetErrorMeanAndImage(bestMappingSubbedA,hB,wB,hA,wA,br_boundary_to_ignore,rgbA,rgbB,width , 0);
        else
            [GCKrms1,regularError] = GetErrorMeanAndImage(bestMappingSubbedA,hB,wB,hA,wA,br_boundary_to_ignore,A,B,width , 1);
        end
        annErrorsA(tableInd+1) = regularError;
        
        debugShow(insideInfo, 'regularError',regularError);
        
        [gckError,bestErrorsNewA] = GetGCKError(A,B,TLboundary,double(width),classType,bestMappingSubbedA,maxKernels,br_boundary_to_ignore);
        if (useTics)
            innerTimes.CandidateCreation = innerTimes.CandidateCreation + toc;
        end
        
    end % if (run_C_code)
    
end % End of Huge Loop

% debug
debugShow(insideInfo, 'bitCountPerTable',bitCountPerTable);
debugShow(insideInfo, 'numKernelsPerAlternative',numKernelsPerAlternative);

% this calculates all filter images (returned in INT16 format)
if (useTics)
    PostProcessTime = tic;
end

if (~CalcErrorImages) && (k > 1) % If in debug mode, the resahpe for final output is done in debug error calculation
    nBestMapping32uA = reshape(nBestMapping32uA', hA,wA,k);
    if (calcBnn)
        nBestMapping32uB = reshape(nBestMapping32uB', hB,wB,k);
    end
end


%% H] FINALIZE - arrange results in outputs
for i = 1 : k
    if (run_C_code)
        [AnnA2B(:,:,2,i),AnnA2B(:,:,1,i)] = ind2sub([hB wB],nBestMapping32uA(:,:,i)+1);
        if (CalcErrorImages)
            AnnA2B(:,:,3,i) = annErrorsA(end);
        end
        if (calcBnn)
            [AnnB2A(:,:,2,i),AnnB2A(:,:,1,i)] = ind2sub([hA wA],nBestMapping32uB(:,:,i)+1);
            if (CalcErrorImages)
                AnnB2A(:,:,3,i) = annErrorsB(end);
            end
        end
    else
        error('this code is obsolete');
        [AnnA2B(:,:,2),AnnA2B(:,:,1)] = ind2sub([hB wB],bestMappingNewA);
        if (CalcErrorImages)
            AnnB2A(:,:,3) = annErrorsB(end);
        end
    end
end

if (patch_mode)
    AnnA2B = AnnA2B(1:end-width+1,1:end-width+1,:,:); % we don't wont the boundaries in this case
end

if (useTics)
    innerTimes.PostProcessTime = toc(PostProcessTime);
    
    % Calcuate total time
    innerTimes.totalTime = sum([innerTimes.PreProcessTime,innerTimes.GCKTime,innerTimes.MemoryConvertTime,innerTimes.CodesTime,...
        innerTimes.HashingTime,innerTimes.CandidateCreation,innerTimes.PostProcessTime]);
end

return

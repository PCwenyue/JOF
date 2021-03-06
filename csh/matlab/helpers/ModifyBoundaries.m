function [currFiltImgs_A, currFiltImgs_B] = ...
    ModifyBoundaries(currFiltImgs_A, currFiltImgs_B,largeVal,calcBnn,hA,hB,width,descriptor_mode , mask)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LargeValDiv2 = floor(largeVal / 2);
if (descriptor_mode == 0)
    % Loop on first size kernels - no image can return these results for projections -
    if (~calcBnn)
    for p=1:length(currFiltImgs_B)% 6
        if (p > 1)
            largeVal = LargeValDiv2;
        end
        currFiltImgs_B{p}(1:hB) = largeVal; % intmax('int16');
        currFiltImgs_B{p}(:, end-width+1:end) = largeVal; % intmax('int16');
        currFiltImgs_B{p}(end-width+1:end, :) = largeVal; % intmax('int16');
        if (~isempty(mask))
            currFiltImgs_B{p}(mask == 1) = largeVal;
        end
    end
    % same for B
    else
        for p=1:length(currFiltImgs_A)% 6
        if (p > 1)
            largeVal = LargeValDiv2;
        end
            currFiltImgs_A{p}(1:hA) = largeVal; % intmax('int16');
            currFiltImgs_A{p}(:, end-width+1:end) = largeVal; % intmax('int16');
            currFiltImgs_A{p}(end-width+1:end, :) = largeVal; % intmax('int16');
            if (~isempty(mask))
                currFiltImgs_A{p}(mask == 1) = largeVal;
            end
        end
    end
else
    if (~calcBnn)
    currFiltImgs_B(1:hB,1,1) = largeVal;
    currFiltImgs_B(:, end-width+1:end,1) = largeVal;
    currFiltImgs_B(end-width+1:end, :,1) = largeVal;
    if (~isempty(mask))
        currFiltImgs_B(mask == 1,1) = largeVal;
    end
    currFiltImgs_B(1:hB,1,2:end) = largeVal;
    currFiltImgs_B(:, end-width+1:end,2:end) = largeVal;
    currFiltImgs_B(end-width+1:end, :,2:end) = largeVal;
    if (~isempty(mask))
        currFiltImgs_B(mask == 1,2:end) = largeVal;
    end
    % same for B
    else
        currFiltImgs_A(1:hA,1,1) = largeVal;
        currFiltImgs_A(:, end-width+1:end,1) = largeVal;
        currFiltImgs_A(end-width+1:end, :,1) = largeVal;
        if (~isempty(mask))
            currFiltImgs_A(mask == 1,1) = largeVal;
        end
        currFiltImgs_A(1:hA,1,2:end) = largeVal;
        currFiltImgs_A(:, end-width+1:end,2:end) = largeVal;
        currFiltImgs_A(end-width+1:end, :,2:end) = largeVal;
        if (~isempty(mask))
            currFiltImgs_A(mask == 1,2:end) = largeVal;
        end
    end
end

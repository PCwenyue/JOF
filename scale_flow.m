function maxrad = scale_flow(flo,maxflow)
% scale flow for later compuation
    maxradu = abs(flo(:,:,1));
    if exist('maxflow','var')
        maxradu(maxradu>maxflow) = 0;
    end
    maxradu = max(max(maxradu));
     
    maxradv = abs(flo(:,:,2));
    if exist('maxflow','var')
        maxradv(maxradv>maxflow) = 0;
    end
    maxradv = max(max(maxradv));
    
    maxrad = sqrt(maxradu^2+maxradv^2);
end
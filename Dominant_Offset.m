function [histmat top_uv] = Dominant_Offset(uv,valid_flag)

[m,n,c] = size(uv);
K = 40; %log_flag = 1;

u = uv(:,:,1); u=u(:)';
v = uv(:,:,2); v=v(:)';

if exist('valid_flag','var') && ~isempty(valid_flag)
    ind = find(valid_flag==1);
    u=u(ind); v=v(ind);
end

maxvalue = max( max(abs(u)), max(abs(v)) );
ui = linspace(-maxvalue,maxvalue,2*maxvalue+1);
vi = linspace(-maxvalue,maxvalue,2*maxvalue+1);

histmat  = hist2(u, v, ui, vi);
histmat = histmat';
% 
% if log_flag == 1
%     histmat=log(histmat+1);
% end

[U,V] = meshgrid(ui,vi);
%figure,surf(U,V,(histmat).^0.5);

top_uv = dominant_uv(histmat,U,V,K);

end

function uv = dominant_uv(histmat,U,V,K)
if ~exist('K','var')
    K = 15;
end

    histmat = histmat(:);
    idx = find(histmat~=0);
    K = min(K,length(idx));
    
    U = U(:)';
    V = V(:)';
    
    [b,ix] = sort(histmat,'descend');
    U = U(ix);
    V = V(ix);
    uv = [U;V;histmat(ix)'];
    uv = uv(:,1:K);
end
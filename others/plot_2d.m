function plot_2d(uv,maxflow,elim_flag)
if ~exist('maxflow','var')
    maxflow=20;
end
if ~exist('elim_flag','var')
    elim_flag = 1;
end

log_flag = 1;

figure,
u = uv(:,:,1); u=u(:)';
v = uv(:,:,2); v=v(:)';

if elim_flag == 1
    ind = find(u<maxflow);
    u=u(ind);
    v=v(ind);
end

n = maxflow*2+1;
ui = linspace(-maxflow,maxflow,n);
vi = linspace(-maxflow,maxflow,n);

ur = interp1(ui,1:numel(ui),u,'nearest');
vr = interp1(vi,1:numel(vi),v,'nearest');

if log_flag
    ur = round(log(ur+1));
    vr = round(log(vr+1));
end

z = accumarray([ur' vr'],1,[n n]);

surf(z);

end
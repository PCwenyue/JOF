function demo_flow(top_homo,m,n)

[x,y] = meshgrid(1:n,1:m);

k = floor(sqrt(length(top_homo)))+1;
figure,
for i=1:length(top_homo)
    tmp_homo = top_homo(i).matrix;
    x2 = tmp_homo(1,1)*x + tmp_homo(1,2)*y + tmp_homo(1,3);
    y2 = tmp_homo(2,1)*x + tmp_homo(2,2)*y + tmp_homo(2,3);
    z2 = tmp_homo(3,1)*x + tmp_homo(3,2)*y + tmp_homo(3,3);
    x2 = x2./z2; y2 = y2./z2;
    
    u = x2 - x;    v = y2 - y;
    uv = cat(3, u, v);
    
    testimg1 = flowToColor(uv); %imshow(testimg)
    subplot(k, k, i),imshow(testimg1)
end

end
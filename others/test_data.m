% funtion to test data term

[error_map,I] = min(error_map,[],1);
I = reshape(I,[m n]);

u = zeros(m,n);
v = zeros(m,n);

for i=1:Knn_count
    tmp_l = find(I==i);
    u(tmp_l) = dominant_uv(1,i);
    v(tmp_l) = dominant_uv(2,i);
end
uv = cat(3,u,v);

testimg1 = flowToColor(uv); %imshow(testimg)
figure,imshow(testimg1)
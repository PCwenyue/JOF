function [pca_im1,pca_im2]=pca_img(im1,im2,patch_wid)

[m,n,c] = size(im1);
pca_dim = max(round(patch_wid*patch_wid*c/7),10);

%% collect data
im1_col = zeros(patch_wid*patch_wid*c,(m-patch_wid+1)*(n-patch_wid+1)); % k x N
for k=1:c
    im1_col( (patch_wid*patch_wid)*(k-1)+1 : patch_wid*patch_wid*k , : ) =...
        im2col(im1(:,:,k),[patch_wid patch_wid],'sliding');
end

im2_col = zeros(patch_wid*patch_wid*c,(m-patch_wid+1)*(n-patch_wid+1)); % k x N
for k=1:c
    im2_col( (patch_wid*patch_wid)*(k-1)+1 : patch_wid*patch_wid*k , : ) =...
        im2col(im2(:,:,k),[patch_wid patch_wid],'sliding');
end

%% PCA
pca_data_num = min(5000,size(im1_col,2));
m1 = mean(im1_col,2);
m2 = mean(im2_col,2);
mean_vec = (m1+m2)/2;
rand_ind1 = randperm( (m-patch_wid+1)*(n-patch_wid+1) );
data1 = im1_col(:,rand_ind1(1:pca_data_num));
rand_ind2 = randperm( (m-patch_wid+1)*(n-patch_wid+1) );
data2 = im2_col(:,rand_ind2(1:pca_data_num));
pca_data = [data1 data2] - repmat(mean_vec,[1 pca_data_num*2]);
cov_matrix = pca_data*pca_data';
[u,s,v] = svd(cov_matrix);
sub_space = v(:,1:pca_dim)';
im1_col = sub_space * im1_col;
im2_col = sub_space * im2_col;

%% PCA image
pca_im1 = zeros(m-patch_wid+1,n-patch_wid+1,pca_dim);
pca_im2 = zeros(m-patch_wid+1,n-patch_wid+1,pca_dim);

for i=1:pca_dim
    tmp = reshape(im1_col(i,:),[m-patch_wid+1 n-patch_wid+1]);
    pca_im1(:,:,i) = tmp;
    tmp = reshape(im2_col(i,:),[m-patch_wid+1 n-patch_wid+1]);
    pca_im2(:,:,i) = tmp;
end

end
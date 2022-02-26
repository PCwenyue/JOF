function top_homo = dense_match(uv,shift)
if~exist('shift','var')
    shift = 0;
end
[m,n,c]=size(uv);
[X1,Y1] = meshgrid(1+shift:n+shift,1+shift:m+shift);

X2 = X1 + uv(:,:,1);
Y2 = Y1 + uv(:,:,2);

% parameter setting
rand_num = 250;thr_pt = 3000;

% sift matches between images
% X1 = pt_match(:,1:2)';
% X2 = pt_match(:,3:4)';
% X1(3,:) = 1;

% homography computing
A = zeros(8,8);
b = zeros(8,1);
tmp_homo_label = 1;
tmp_pt_match = zeros(4,4);
for x_range = 1:4
    for y_range = 1:4
        [x_range y_range]
        for i = 1:rand_num
            x_rand = randperm(floor(n/4)); x_rand = x_rand + (x_range-1)*floor(n/4);    x_rand = x_rand(1:4);
            y_rand = randperm(floor(m/4)); y_rand = y_rand + (y_range-1)*floor(m/4);    y_rand = y_rand(1:4);
            for j=1:4
                tmp_pt_match(j,1) = X1(y_rand(j),x_rand(j)); %x1
                tmp_pt_match(j,2) = Y1(y_rand(j),x_rand(j)); %y1
                tmp_pt_match(j,3) = X2(y_rand(j),x_rand(j)); %x2
                tmp_pt_match(j,4) = Y2(y_rand(j),x_rand(j)); %y2
            end

            % computing homography
            A(:) = 0; b(:) = 0;
            for j = 1:4
                x1 = tmp_pt_match(j,1);        y1 = tmp_pt_match(j,2);
                x2 = tmp_pt_match(j,3);        y2 = tmp_pt_match(j,4);

                tmpA = [x1 y1 1 0  0  0 -x2*x1 -x2*y1;...
                        0  0  0 x1 y1 1 -y2*x1 -y2*y1];
                tmpb = [x2; y2];
                A(j*2-1:j*2,:) = tmpA; b(j*2-1:j*2) = tmpb;
            end
            tmphomo = A\b;
            tmphomo = [tmphomo(1) tmphomo(2) tmphomo(3);...
                       tmphomo(4) tmphomo(5) tmphomo(6);...
                       tmphomo(7) tmphomo(8) 1];

            % computing residual
            newX2 = tmphomo(1,1) * X1 + tmphomo(1,2) * Y1 + tmphomo(1,3);
            newY2 = tmphomo(2,1) * X1 + tmphomo(2,2) * Y1 + tmphomo(2,3);
            newZ2 = tmphomo(3,1) * X1 + tmphomo(3,2) * Y1 + tmphomo(3,3);
            newX2 = newX2./newZ2;    newY2 = newY2./newZ2;

            tmp_res_X = abs(X2 - newX2); tmp_res_X=(tmp_res_X<1);
            tmp_res_Y = abs(Y2 - newY2); tmp_res_Y=(tmp_res_Y<1);
            tmp_res = tmp_res_X .* tmp_res_Y;

            tmp_valid = find(tmp_res);

            if length(tmp_valid) > thr_pt
                if tmp_homo_label > 1 % check for homography repetition
                    repeat_flag = 0;
                    for k=1:tmp_homo_label-1
                        label1 = top_homo(k).label;
                        label2 = tmp_valid;
                        if length(intersect(label1,label2)) > 0.6 *min(length(label1),length(label2));
                            repeat_flag = 1;
                            break;
                        end
                    end
                end
                if tmp_homo_label == 1 || repeat_flag ==0
                    top_homo(tmp_homo_label).matrix = tmphomo;
                    top_homo(tmp_homo_label).label = tmp_valid;
                    tmp_homo_label = tmp_homo_label+1;
                else % repetition occurs
                    if length(tmp_valid) > length(top_homo(k).label)
                        top_homo(k).matrix = tmphomo;
                        top_homo(k).label = tmp_valid;
                    end
                end
            end
        end
    end
end

end
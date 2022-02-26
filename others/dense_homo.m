function top_homo = dense_homo(uv,shift,valid_flag)
% give every coordinate a correction?
if~exist('shift','var')||isempty(shift)
    shift = 0;
end
[m,n,c]=size(uv);
if~exist('valid_flag','var')||isempty(valid_flag)
    valid_flag = ones(m,n);
end

[X1,Y1] = meshgrid(1+shift:n+shift,1+shift:m+shift);

X2 = X1 + uv(:,:,1);Y2 = Y1 + uv(:,:,2);

% parameter setting
rand_num = 250;
thr_pt = 3000;
pt_ransac = 3;
pyramid = 2;

% homography computing
A = zeros(pt_ransac*2,4);b = zeros(pt_ransac*2,1);
tmp_homo_label = 1;
tmp_pt_match = zeros(pt_ransac,4);

for x_range = 1:pyramid
    for y_range = 1:pyramid
        %% for each pyramid
        [x_range y_range]
        for i = 1:rand_num % how many ransac iterations for each pyramid
            %i
            tmp_valid_flag = valid_flag( (y_range-1)*floor(m/pyramid)+1 : y_range*floor(m/pyramid) ,...
                (x_range-1)*floor(n/pyramid)+1 : x_range*floor(n/pyramid) );
            valid_index = find(tmp_valid_flag);
            
            if length(valid_index) > 2000
                % tmp local
                tmp_X1 = X1( (y_range-1)*floor(m/pyramid)+1 : y_range*floor(m/pyramid) ,...
                (x_range-1)*floor(n/pyramid)+1 : x_range*floor(n/pyramid) );
                tmp_X2 = X2( (y_range-1)*floor(m/pyramid)+1 : y_range*floor(m/pyramid) ,...
                (x_range-1)*floor(n/pyramid)+1 : x_range*floor(n/pyramid) );
                tmp_Y1 = Y1( (y_range-1)*floor(m/pyramid)+1 : y_range*floor(m/pyramid) ,...
                (x_range-1)*floor(n/pyramid)+1 : x_range*floor(n/pyramid) );
                tmp_Y2 = Y2( (y_range-1)*floor(m/pyramid)+1 : y_range*floor(m/pyramid) ,...
                (x_range-1)*floor(n/pyramid)+1 : x_range*floor(n/pyramid) );
                
                tmprand = randperm(length(valid_index));
                valid_index = valid_index(tmprand(1:pt_ransac));
                %x_rand = randperm(floor(n/pyramid)); x_rand = x_rand + (x_range-1)*floor(n/pyramid);    %x_rand = x_rand(1:pt_ransac);
                %y_rand = randperm(floor(m/pyramid)); y_rand = y_rand + (y_range-1)*floor(m/pyramid);    %y_rand = y_rand(1:pt_ransac);
                %while count < 4
                %    if valid_flag(y_rand(j),x_rand(j))==1
                for count=1:pt_ransac
                    tmp_pt_match(count,1) = tmp_X1(valid_index(count)); %x1
                    tmp_pt_match(count,2) = tmp_Y1(valid_index(count)); %y1
                    tmp_pt_match(count,3) = tmp_X2(valid_index(count)); %x2
                    tmp_pt_match(count,4) = tmp_Y2(valid_index(count)); %y2
                end

                % counting the triangle area
                aa = sqrt((tmp_pt_match(1,1)-tmp_pt_match(2,1))^2+(tmp_pt_match(1,2)-tmp_pt_match(2,2))^2);%d12
                bb = sqrt((tmp_pt_match(2,1)-tmp_pt_match(3,1))^2+(tmp_pt_match(2,2)-tmp_pt_match(3,2))^2);%d23
                cc = sqrt((tmp_pt_match(3,1)-tmp_pt_match(1,1))^2+(tmp_pt_match(3,2)-tmp_pt_match(1,2))^2);%d13
                ss = (aa+bb+cc)/2;
                pix_num = sqrt(ss*(ss-aa)*(ss-bb)*(ss-cc));

                % computing homography
                A(:) = 0; b(:) = 0;
                for j = 1:pt_ransac
                    x1 = tmp_pt_match(j,1);y1 = tmp_pt_match(j,2);
                    x2 = tmp_pt_match(j,3);y2 = tmp_pt_match(j,4);

                    tmpA = [x1 -y1 1 0;...
                            y1 x1 0 1];
                    tmpb = [x2; y2];
                    A(j*2-1:j*2,:) = tmpA; b(j*2-1:j*2) = tmpb;
                end
                tmphomo = A\b;
                tmphomo = [tmphomo(1) -tmphomo(2) tmphomo(3);...
                           tmphomo(2) tmphomo(1) tmphomo(4)];

                % computing residual
                newX2 = tmphomo(1,1) * X1 + tmphomo(1,2) * Y1 + tmphomo(1,3);
                newY2 = tmphomo(2,1) * X1 + tmphomo(2,2) * Y1 + tmphomo(2,3);

                tmp_res_X = abs(X2 - newX2); tmp_res_X=(tmp_res_X<.5);
                tmp_res_Y = abs(Y2 - newY2); tmp_res_Y=(tmp_res_Y<.5);
                tmp_res = tmp_res_X .* tmp_res_Y;
                tmp_res = tmp_res .* valid_flag;
                tmp_valid = find(tmp_res);

                if length(tmp_valid) > max(thr_pt,pix_num)
                    tx1 = X1(tmp_valid); ty1 = Y1(tmp_valid);
                    tx2 = X2(tmp_valid); ty2 = Y2(tmp_valid);
                    tmphomo2 = homo_refinement(tx1,ty1,tx2,ty2);
                    top_homo(tmp_homo_label).matrix = tmphomo2;
                    top_homo(tmp_homo_label).label = length(tmp_valid);
                    tmp_homo_label = tmp_homo_label + 1;
                    valid_flag(tmp_valid) = 0; % without replacement
                end
            else
                break;
            end
        end
    end
end

if ~exist('top_homo','var')
    top_homo=[];
end

end

%% homography refinement
function tmphomo=homo_refinement(x1,y1,x2,y2)
    count = length(x1);
    A = zeros(2*count,4); b = zeros(2*count,1);
    for j = 1:count
        tmpA = [x1(j) -y1(j) 1 0;...
                y1(j)  x1(j) 0 1];
        tmpb = [x2(j); y2(j)];
        A(j*2-1:j*2,:) = tmpA; b(j*2-1:j*2) = tmpb;
    end
    tmphomo = A\b;
    tmphomo = [tmphomo(1) -tmphomo(2) tmphomo(3);...
               tmphomo(2) tmphomo(1) tmphomo(4)];
end

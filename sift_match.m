function top_homo = sift_match(I0,I1)

% parameter setting
rand_num = 1000;
thr_pt = 100;

% sift matches between images
pt_match = siftmatch(I0, I1); % [Nx4], x1,y1,x2,y2
delete('tmp.key');    delete('tmp.pgm');
match_num = size(pt_match,1);

X1 = pt_match(:,1:2)';
X2 = pt_match(:,3:4)';
X1(3,:) = 1;

% homography computing by RANSAC
A = zeros(8,8);
b = zeros(8,1);
tmp_homo_label = 1;
for i = 1:rand_num
    a = randperm(match_num);
    a = a(1:4);
    
    tmp_pt_match = pt_match(a,:);
    
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
    tmpX2 = tmphomo * X1;
    tmpX2 = tmpX2 ./ repmat(tmpX2(3,:),[3 1]);
    tmp_res = X2 - tmpX2(1:2,:);
    
    tmp_valid = max(abs(tmp_res),[],1);
    tmp_valid = find(tmp_valid<1);
    
    if length(tmp_valid) > thr_pt
        if tmp_homo_label > 1 % check for homography repetition
            repeat_flag = 0;
            for k=1:tmp_homo_label-1
                label1 = top_homo(k).label;
                label2 = tmp_valid;
                if length(intersect(label1,label2)) > 0.6 * min(length(label1),length(label2));
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

%% homography refinement
if exist('top_homo','var')
    for i = 1:length(top_homo)
        tmp_pt = pt_match(top_homo(i).label,:);%x1,y1,x2,y2
        tmp_num = size(tmp_pt,1);
        A = zeros(2*tmp_num,8);    b = zeros(2*tmp_num);
        for j = 1:tmp_num
            x1 = tmp_pt(j,1);        y1 = tmp_pt(j,2);
            x2 = tmp_pt(j,3);        y2 = tmp_pt(j,4);

            tmpA = [x1 y1 1 0  0  0 -x2*x1 -x2*y1;...
                    0  0  0 x1 y1 1 -y2*x1 -y2*y1];
            tmpb = [x2; y2];
            A(j*2-1:j*2,:) = tmpA; b(j*2-1:j*2) = tmpb;
        end
        tmphomo = A\b;
        tmphomo = [tmphomo(1) tmphomo(2) tmphomo(3);...
                   tmphomo(4) tmphomo(5) tmphomo(6);...
                   tmphomo(7) tmphomo(8) 1];
        top_homo(i).matrix = tmphomo;
    end
else
    top_homo=[];
end

end
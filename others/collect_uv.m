function uv_core=collect_uv(uv)

uv = uv(1:2,:);

% uv_core = round(uv/3);
% uv_core = uv_core(1,:) + uv_core(2,:)/100;
% uv_core = unique(uv_core);
% 
% uv_core(2,:) = (uv_core-round(uv_core))*100;
% uv_
uv_core = [];
k = 1;
for i=1:size(uv,2)
    tmp_uv = uv(:,i);
    
    repeat_flag = 0;
    for j=1:k-1
        tmp_uv_core=uv_core(j).grid;
        if max(abs(tmp_uv_core-tmp_uv)) <2
            repeat_flag = 1;
            uv_core(j).uv = [uv_core(j).uv tmp_uv];
            break;
        end
    end
    
    if repeat_flag == 0
        uv_core(k).grid = round(tmp_uv/3)*3;
        uv_core(k).uv = tmp_uv;
        k = k+1;
    end
end

end
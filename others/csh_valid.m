function flag = csh_valid(uv,top_uv)
% only mark top velocities in uv as 1 in flag

%% checking only the top_velocity component
[m,n,c] = size(uv);
flag = zeros(m,n);

for i=1:size(top_uv,2)
    tmp_u = top_uv(1,i);
    tmp_v = top_uv(2,i);
    tmp_flag_u = (uv(:,:,1)==tmp_u);
    tmp_flag_v = (uv(:,:,2)==tmp_v);
    
    tmp_flag = tmp_flag_u .* tmp_flag_v;
    flag = flag + tmp_flag;
end

end
function top_homo3 = homo_pick(top_homo2,homo_count)

count = length(top_homo2);
homo_num = zeros(count,1);

for i=1:count
    if length(top_homo2(i).label)>1
        homo_num(i) = length(top_homo2(i).label);
    else
        homo_num(i) = top_homo2(i).label;
    end
end
[b,ix] = sort(homo_num,'descend');
ix = ix(1:homo_count);

top_homo3 = top_homo2(ix);

end
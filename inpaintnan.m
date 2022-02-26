function final_img=inpaintnan(img)

inpaint_wid = 5;
a = img;a(~isnan(a)) = 1; a(isnan(a)) = 0;

%[m,n]=size(img);
%b = edge(a); b = find(b);
%bx = floor((b-1)/m)+1; by = mod(b-1,m)+1;

img(isnan(img))=0;
final_img = img;
final_img(isnan(final_img))=0;

% moving left until no update available
update_flag = 1;
while update_flag == 1
    new_img = [final_img(:,inpaint_wid+1:end,:) final_img(:,end-inpaint_wid+1:end,:)];
    new_a = [a(:,inpaint_wid+1:end,:) a(:,end-inpaint_wid+1:end,:)];
    new_a = new_a - a;    new_a(new_a~=1) = 0;
    if sum(sum(new_a))==0
        break;
    else
        a = a + new_a;
        final_img = final_img + new_a.*new_img;
    end
end

% moving right until no update available
update_flag = 1;
while update_flag == 1
    new_img = [final_img(:,1:inpaint_wid,:) final_img(:,1:end-inpaint_wid,:)];
    new_a = [a(:,1:inpaint_wid,:) a(:,1:end-inpaint_wid,:)];
    new_a = new_a - a;    new_a(new_a~=1) = 0;
    if sum(sum(new_a))==0
        break;
    else
        a = a + new_a;
        final_img = final_img + new_a.*new_img;
    end
end

% moving up until no update available
update_flag = 1;
while update_flag == 1
    new_img = [final_img(inpaint_wid+1:end,:,:); final_img(end-inpaint_wid+1:end,:,:)];
    new_a = [a(inpaint_wid+1:end,:,:); a(end-inpaint_wid+1:end,:,:)];
    new_a = new_a - a;    new_a(new_a~=1) = 0;
    if sum(sum(new_a))==0
        break;
    else
        a = a + new_a;
        final_img = final_img + new_a.*new_img;
    end
end

% moving down until no update available
update_flag = 1;
while update_flag == 1
    new_img = [final_img(1:inpaint_wid,:,:); final_img(1:end-inpaint_wid,:,:)];
    new_a = [a(1:inpaint_wid,:,:); a(1:end-inpaint_wid,:,:)];
    new_a = new_a - a;    new_a(new_a~=1) = 0;
    if sum(sum(new_a))==0
        break;
    else
        a = a + new_a;
        final_img = final_img + new_a.*new_img;
    end
end

end
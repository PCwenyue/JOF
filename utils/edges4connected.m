function E = edges4connected(height,width,pair_flag)

% EDGE4CONNECTED Creates edges where each node
%   is connected to its four adjacent neighbors on a 
%   height x width grid.
%   E - a vector in which each row i represents an edge
%   E(i,1) --> E(i,2). The edges are listed is in the following 
%   neighbor order: down,up,right,left, where nodes 
%   indices are taken column-major.
%
%   (c) 2008 Michael Rubinstein, WDI R&D and IDC
%   $Revision$
%   $Date$
%
if ~exist('pair_flag','var')
    pair_flag = 0;
end
N = height*width;
I = []; J = [];
% connect vertically (down, then up)
is = [1:N]'; is([height:height:N])=[];
js = is+1;
if pair_flag == 0
    I = [I;is];
    J = [J;js];
else
    I = [I;is;js];
    J = [J;js;is];
end
% connect horizontally (right, then left)
is = [1:N-height]';
js = is+height;
if pair_flag == 0
    I = [I;is];
    J = [J;js];
else
    I = [I;is;js];
    J = [J;js;is];
end
E = [I,J];

end
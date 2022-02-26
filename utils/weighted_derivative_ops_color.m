function [weights_x weights_y] = weighted_derivative_ops_color(M, N, I0, beta, nu, zero_pedding)

if ~exist('zero_pedding','var')
    zero_pedding = 1;
end
	%% Convert RGB image to LAB
	%I0LAB = vl_xyz2lab(vl_rgb2xyz(double(I0)/255)); 
    I0LAB = RGB2Lab(double(I0/255));
    I0RGB = double(I0);

	MN = M*N;
	[X,Y]  = meshgrid(1:N, 1:M);

	%% Derivative in X direction
	edges = [ vec(sub2ind([M N], Y(:, 1:end-1), X(:, 1:end-1))) ...
	          vec(sub2ind([M N], Y(:, 1:end-1), X(:, 1:end-1)+1)) ];
	nedges = size(edges, 1);

	%% Compute image features
	weights = makeweights(edges, [vec(I0LAB(:,:,1)) vec(I0LAB(:,:,2)) vec(I0LAB(:,:,3))], beta);

	weights_x = nu + (1-nu)*weights;
    if zero_pedding == 1
        weights_x = [ zeros(M,1) ; weights_x ];
    end
	%% The derivative operator in X direction
	%D1 = sparse([edges(:,1); edges(:, 1)], [edges(:, 1); edges(:,2)], [-weights_; weights_], MN, MN);

	%% Derivative in Y direction
	edges = [vec(sub2ind([M N], Y(1:end-1,:), X(1:end-1,:)))  ...
	         vec(sub2ind([M N], Y(1:end-1,:)+1, X(1:end-1,:)))];
					
	nedges = size(edges, 1);

	%% Compute image features
	weights = makeweights(edges, [vec(I0LAB(:,:,1)) vec(I0LAB(:,:,2)) vec(I0LAB(:,:,3))], beta);

	weights_y = nu + (1-nu)*weights;  
    
    if zero_pedding == 1
        weights_yy = zeros(MN,1);
        for i=1:N
            weights_yy( (i-1)*M+2:i*M ) = weights_y( (i-1)*(M-1)+1:i*(M-1) );
        end
        weights_y = weights_yy;
    end

	%% The derivative operator in Y direction
	%D2 = sparse([edges(:,1); edges(:, 1)], [edges(:, 1); edges(:,2)], [-weights_; weights_], MN, MN);

end
function f = vec(f)
f = f(:);
end
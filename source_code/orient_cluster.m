function theta = orient_cluster(I, J, index_cluster, c_x, c_y, height)

% size_fig = size(fig_cluster);
% fig = zeros(size_fig);
theta = 0;

[r c] = size(I);
X = zeros(2,r);

X(1,:) = (height-I) - c_y;
X(2,:) = J - c_x;

cov_mat = X*X';
[V] = eig(cov_mat);

y = V(1); %+c_y;
%y = size_fig(1)-y;
x = V(2); %+c_x;

% theta = atan(y/x);
% theta = pi - theta;

theta = atan(x/y);
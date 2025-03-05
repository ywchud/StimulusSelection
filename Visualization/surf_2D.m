function surf_2D(X,Y,Z)
%Darren
%replaces surf() to plot a 2D map with axis X and Y of matrix Z assigned to a colorbar
%An extra row and column is added to the matrix Z to make the map
%size identical to the matrix
%similarly to surf, X and Y should have the same size as Z, check meshgrid()
%X and Y should contain the edges of each bin, excluding the last index e.g. (1:end-1)
Z=[Z;Z(end,:)];
Z=[Z Z(:,end)];
dY=median(diff(Y(:,1)));
dX=median(diff(X(1,:)));
Y=[Y;Y(end,:)+dY];
Y=[Y Y(:,end)];
X=[X;X(end,:)];
X=[X X(:,end)+dX];
surf(X,Y,Z)
view(0,90)
end
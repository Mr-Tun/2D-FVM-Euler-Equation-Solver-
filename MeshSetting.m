function [Xc,Yc,omega,SI,SJ,right,upper,left,lower] = MeshSetting(X,Y)
%MESHSETTING is used to get the grid parameter by the coordinate of grid
%points for the FVM
[nx, ny] = size(X);
omega = zeros(nx-1,ny-1);
SI = zeros(nx,ny-1);
SJ = zeros(nx-1,ny);
right = zeros(nx-1,ny-1,2);
upper = zeros(nx-1,ny-1,2);
left = zeros(nx-1,ny-1,2);
lower = zeros(nx-1,ny-1,2);
for i=1:nx-1
    for j=1:ny-1
        %        4----3
        %       /    /
        %      /    /
        %     1----2
        % use 4 points of the cell to compute the Omega, S, n vector
        x1 = X(i, j);
        y1 = Y(i, j);
        x2 = X(i+1, j);
        y2 = Y(i+1, j);
        x3 = X(i+1, j+1);
        y3 = Y(i+1, j+1);
        x4 = X(i, j+1);
        y4 = Y(i, j+1);
        % the length of the cell I,J
        SI(i, j) = sqrt((x1 - x4)^2 + (y1 - y4)^2);
        SJ(i, j) = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        % the normal identity vector of the cell I,J 

        upper(i, j,:) = -[y3-y4,-x3+x4]/norm([y3-y4,-x3+x4]);
        right(i, j,:) = -[y2-y3,-x2+x3]/norm([y2-y3,-x2+x3]);
        lower(i, j,:) = -[y1-y2,-x1+x2]/norm([y1-y2,-x1+x2]);
        left(i, j,:) =  -[y4-y1,-x4+x1]/norm([y4-y1,-x4+x1]); 
%         upper(i, j,:) = -[y3-y4,-x3+x4];
%         right(i, j,:) = -[y2-y3,-x2+x3];
%         lower(i, j,:) = -[y1-y2,-x1+x2];
%         left(i, j,:) =  -[y4-y1,-x4+x1]; 
        if j==ny-1
            SJ(i,j+1) = sqrt((x3 - x4)^2 + (y3 - y4)^2);
        end
        if i==nx-1
            SI(i+1,j) = sqrt((x3 - x2)^2 + (y3 - y2)^2);
        end
        % the area of triangle 1: point 1,2,3
        omega1 = 0.5 * abs((x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1));
        % the area of triangle 2: 1,3,4
        omega2 = 0.5 * abs((x3 - x1) * (y4 - y1) - (x4 - x1) * (y3 - y1));
        % the Area of the quadrilateral grid
        omega(i, j) = omega1 + omega2;
    end
end
 Xc = (X(1:end-1, 1:end-1) + X(2:end, 1:end-1) + X(1:end-1, 2:end) + X(2:end, 2:end)) / 4;
 Yc = (Y(1:end-1, 1:end-1) + Y(2:end, 1:end-1) + Y(1:end-1, 2:end) + Y(2:end, 2:end)) / 4;
end


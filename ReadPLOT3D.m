function [X,Y]=ReadPLOT3D(filename)
%% open the file
fid = fopen(filename, 'r');
    if fid == -1
        error('can open the file: %s', filename);
    end
%% read Imax and Jmax of grid points
    dims = fscanf(fid, '%d', 2);
    Imax = dims(1);
    Jmax = dims(2);    
%% read the location of the grid points
    X = fscanf(fid, '%f', [Imax, Jmax]);
    Y = fscanf(fid, '%f', [Imax, Jmax]);
%% close the file
    fclose(fid);
%% depict the grid
X = X*10^(-3); Y = Y*10^(-3);
figure(1)
hold on;
for i = 1:size(X, 1)
    line(X(i, :), Y(i, :), 'Color', 'b'); % 绘制每一行的网格线
end
for j = 1:size(X, 2)
    line(X(:, j), Y(:, j), 'Color', 'b'); % 绘制每一列的网格线
end
hold off;
title('2D Grid of blunt-body');
xlabel('\itx\rm/m');
ylabel('\ity\rm/m');

% Set the scale of x-axis and y-axis same 
axis equal;
grid on;
disp(['the ',filename,' is read']);
end


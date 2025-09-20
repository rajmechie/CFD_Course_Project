% Load the data from the file
data = load('stream_1000.dat'); % Adjust the file path if needed

% Extract x, y, and psi
x = data(:, 1); % First column for x
y = data(:, 2); % Second column for y
psi = data(:, 3); % Third column for psi

% Determine the grid of x and y values
unique_x = unique(x);
unique_y = unique(y);
[X, Y] = meshgrid(unique_x, unique_y);

% Reshape psi to fit the grid
Psi = griddata(x, y, psi, X, Y);

% Plot the streamline pattern using contour
figure;
contour(X, Y, Psi, 500); % 500 contour levels
colorbar; % Optional: add a colorbar to visualize psi values
xlabel('x');
ylabel('y');
title('Vorticity Contour at Re=1000');
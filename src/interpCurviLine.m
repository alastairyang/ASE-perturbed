function [line_equal] = interpCurviLine(line, n_points)
% INTERPCURVILINE Interpolates a curvilinear line to have equally spaced points
%
% Inputs:
%   line     - Nx2 array of [x,y] coordinates defining the line
%   n_points - Number of desired points in the output line
%
% Output:
%   line_equal - nx2 array of equally spaced [x,y] coordinates
%
% Example:
%   line_equal = interpCurviLine(mid_flowline, 1000);

% Input checking
if size(line,2) ~= 2
    error('Input line must be Nx2 array of [x,y] coordinates');
end
if ~isscalar(n_points) || n_points < 2
    error('n_points must be a scalar > 1');
end

% Calculate cumulative distance along the curve
dx = diff(line(:,1));
dy = diff(line(:,2));
segment_lengths = sqrt(dx.^2 + dy.^2);
cum_distance = [0; cumsum(segment_lengths)];

% Create new evenly spaced points along the cumulative distance
new_distances = linspace(0, cum_distance(end), n_points);

% Interpolate x and y coordinates separately using cumulative distance
x_new = interp1(cum_distance, line(:,1), new_distances);
y_new = interp1(cum_distance, line(:,2), new_distances);

% Combine into new equally-spaced line
line_equal = [x_new', y_new'];

end

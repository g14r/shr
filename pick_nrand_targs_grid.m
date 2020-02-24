function [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs_grid(n_targs, varargin)
%% function [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs_grid(n_targs, varargin)
%
% Pick n successive spatial targets randomly positioned within the range
% defined by x_range ([x_min x_max]) and y_range ([y_min y_max]) for x, y
% coordinates. Return x and y coordinates of each target separately. min_d
% and max_d are the additional optional input constraints that dictate that
% successive targets cannot be closer than a minimum distance, or farther
% than a maximun distance (in % of largest range). The purpose is to
% prevent overlapping targets and ensure that all sequences have roughly
% the same path length.
%
% INPUTS
%     n_targs:    number of targets to be created
%
% VARARGIN
%     x_range:    spatial range of x coordinates [x_min x_max]
%     y_range:    spatial range of y coordinates [y_min y_max]
%     grid_size:  how many lines are defining the grid
%
% OUTPUTS
%     x_coord:    [1 x n_targs] vector of x coordinates for each target
%     y_coord:    [1 x n_targs] vector of y coordinates for each target
%     dist:       [1 x n_targs] vector of distances between successive targets (dist(1) will be 0)
%     tgt_num:    [1 : n_targs] vector of target numbers for each sequence
%
% USAGE
%     [x_coord, y_coord, dist] = pick_nrand_targs_grid(14);
%     [x_coord, y_coord, dist] = pick_nrand_targs_grid(14, 'x_range',[-10 10], 'y_range',[-10 10], 'grid_size',7);
%
% --
% gariani@uwo.ca - 2020.02.18

% deal with eventual variable input arguments and set defaults
x_range = []; y_range = []; grid_size = 7; % how many lines splitting the grid?
vararginoptions(varargin, {'x_range', 'y_range', 'grid_size'});
if isempty(x_range) || isempty(y_range)
    x_range(1,1) = input('What is the minimum x coordinate? ');
    x_range(1,2) = input('What is the maximum x coordinate? ');
    y_range(1,1) = input('What is the minimum y coordinate? ');
    y_range(1,2) = input('What is the maximum y coordinate? ');
end

% create the regular grid of possible target locations
xc = linspace(x_range(1), x_range(2), grid_size);
yc = linspace(y_range(1), y_range(2), grid_size);
x_grid = nan(grid_size^2,1); y_grid = nan(grid_size^2,1);
for xi = 1 : numel(xc)
    for yi = 1 : numel(yc)
        % remove corners
        if (xi==1 && yi==1) ...
                || (xi==grid_size && yi==1) ...
                || (xi==1 && yi==grid_size) ...
                || (xi==grid_size && yi==grid_size)
            continue
        else
            x_grid(numel(xc)*(xi-1) + yi, 1) = xc(xi);
            y_grid(numel(xc)*(xi-1) + yi, 1) = yc(yi);
        end
    end
end
% sanity check:
plt.scatter(x_grid,y_grid, 'regression','none'); axis equal; hold on;
% plot([x(2), x(grid_size*2+1)], [y(2), y(grid_size*2+1)], 'r', 'linewidth',2); hold off;
min_d = round( sqrt( (x_grid(2) - x_grid(grid_size*2+1))^2 + (y_grid(2) - y_grid(grid_size*2+1))^2 ) ,2);
x = unique(x_grid(~isnan(x_grid)));
y = unique(y_grid(~isnan(y_grid)));

% preallocate outputs
x_coord = zeros(1,n_targs);
y_coord = zeros(1,n_targs);
dist    = zeros(1,n_targs);
x_cent  = x_range(1) + diff(x_range)/2; % x of the center;
y_cent  = y_range(1) + diff(y_range)/2; % y of the center
for ct = 1 : n_targs % ct = current target
    % pick random coordinates from the grid options
    x_coord(1,ct) = x(randi(numel(x)));
    y_coord(1,ct) = y(randi(numel(y)));
    if ct ==1
        % start from the center of the range
        x_coord(1,ct) = x_cent;
        y_coord(1,ct) = y_cent;
        while      (x_coord(1,ct)==x(1)   && y_coord(1,ct)==y(1)) ...
                || (x_coord(1,ct)==x(end) && y_coord(1,ct)==y(end)) ...
                || (x_coord(1,ct)==x(1)   && y_coord(1,ct)==y(end)) ...
                || (x_coord(1,ct)==x(end) && y_coord(1,ct)==y(1))
            x_coord(1,ct) = x(randi(numel(x)));
            y_coord(1,ct) = y(randi(numel(y)));
        end
    else
        dist(1,ct) = round( sqrt( (x_coord(1,ct) - x_coord(1,ct-1))^2 + (y_coord(1,ct) - y_coord(1,ct-1))^2 ), 2); % distance between successive targets
        % check if successive targets follow distance constraints
        while      (x_coord(1,ct)==x(1)   && y_coord(1,ct)==y(1)) ...
                || (x_coord(1,ct)==x(end) && y_coord(1,ct)==y(end)) ...
                || (x_coord(1,ct)==x(1)   && y_coord(1,ct)==y(end)) ...
                || (x_coord(1,ct)==x(end) && y_coord(1,ct)==y(1)) ...
                || (dist(1,ct) ~= min_d) ...
                || (x_coord(1,ct)==x_coord(1,ct-1) && y_coord(1,ct)==y_coord(1,ct-1)) ...
                || sum( x_coord(1,ct)==x_coord(1,:) & y_coord(1,ct)==y_coord(1,:) ) > 1 % check whether any target repeats withing the same sequence
            % if any of these conditions are met, keep searching until they are not
            x_coord(1,ct) = x(randi(numel(x)));
            y_coord(1,ct) = y(randi(numel(y)));
            dist(1,ct) = round( sqrt( (x_coord(1,ct) - x_coord(1,ct-1))^2 + (y_coord(1,ct) - y_coord(1,ct-1))^2 ), 2);
        end
    end
end
tgt_num = 1 : n_targs;
end
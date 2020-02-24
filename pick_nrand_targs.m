function [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs(n_targs, varargin)
%% function [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs(n_targs, varargin)
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
%     min_d:      minimum distance between successive targets in percentage of largest range (default is 20%)
%     max_d:      maximum distance between successive targets in percentage of largest range (default is 21%)
%
% OUTPUTS
%     x_coord:    [1 x n_targs] vector of x coordinates for each target
%     y_coord:    [1 x n_targs] vector of y coordinates for each target
%     dist:       [1 x n_targs] vector of distances between successive targets (dist(1) will be 0)
%     tgt_num:    [1 : n_targs] vector of target numbers for each sequence
%
% USAGE
%     [x_coord, y_coord, dist] = pick_nrand_targs(14);
%     [x_coord, y_coord, dist] = pick_nrand_targs(14, 'x_range',[-10 10], 'y_range',[-10 10], 'min_d',20, 'max_d',21);
%
% --
% gariani@uwo.ca - 2020.02.18

% deal with eventual variable input arguments and set defaults
x_range = []; y_range = [];
min_d = []; max_d = [];
vararginoptions(varargin, {'x_range', 'y_range', 'min_d', 'max_d'});
if isempty(x_range) || isempty(y_range)
    x_range(1,1) = input('What is the minimum x coordinate? ');
    x_range(1,2) = input('What is the maximum x coordinate? ');
    y_range(1,1) = input('What is the minimum y coordinate? ');
    y_range(1,2) = input('What is the maximum y coordinate? ');
end
if isempty(min_d) || isempty(max_d)
    min_d = (0.20 * diff(max(x_range, y_range))); % default is 20% of largest range
    max_d = (0.21 * diff(max(x_range, y_range))); % default is 21% of largest range
else
    min_d = ((min_d/100) * diff(max(x_range, y_range))); % expressed in % of largest range
    max_d = ((max_d/100) * diff(max(x_range, y_range))); % expressed in % of largest range
end

% preallocate outputs
x_coord = zeros(1,n_targs);
y_coord = zeros(1,n_targs);
dist    = zeros(1,n_targs);
x_cent  = x_range(1) + diff(x_range)/2; % x of the center;
y_cent  = y_range(1) + diff(y_range)/2; % y of the center
for ct = 1 : n_targs % ct = current target
    % pick random coordinates from a uniform distribution within the range
    x_coord(1,ct) = unifrnd(min(x_range), max(x_range));
    y_coord(1,ct) = unifrnd(min(y_range), max(y_range));
    if ct ==1
        % start from the center of the range
        x_coord(1,ct) = x_cent;
        y_coord(1,ct) = y_cent;
        continue
    else
        dist(1,ct) = sqrt( (x_coord(1,ct) - x_coord(1,ct-1))^2 + (y_coord(1,ct) - y_coord(1,ct-1))^2 ); % distance between successive targets
        dist_all = sqrt( (x_coord(1,ct) - x_coord(1,:)).^2 + (y_coord(1,ct) - y_coord(1,:)).^2 ); % distance of current target from all otehr targets
        dist_cent = sqrt( (x_coord(1,ct) - x_cent)^2 + (y_coord(1,ct) - y_cent)^2 ); % distance of current target from the center (home position)
        while (dist(1,ct) <= min_d || dist(1,ct) >= max_d) ... % check if successive targets follow distance constraints
                || (any(dist_all>0 & dist_all<=(0.05 * diff(max(x_range, y_range))))) ... % check whether distance between any two targets is smaller than 5% of largest range
                || dist_cent <= (0.02 * diff(max(x_range, y_range))) ... % check that no target appears exactly in the home position (center)
                || sum( x_coord(1,ct)==x_coord(1,:) & y_coord(1,ct)==y_coord(1,:) ) > 1 % check whether any target repeats withing the same sequence
            % if any of these conditions are met, keep searching until they are not
            x_coord(1,ct) = unifrnd(min(x_range), max(x_range));
            y_coord(1,ct) = unifrnd(min(y_range), max(y_range));
            dist(1,ct) = sqrt( (x_coord(1,ct) - x_coord(1,ct-1))^2 + (y_coord(1,ct) - y_coord(1,ct-1))^2 );
            dist_all = sqrt( (x_coord(1,ct) - x_coord(1,:)).^2 + (y_coord(1,ct) - y_coord(1,:)).^2 );
            dist_cent = sqrt( (x_coord(1,ct) - x_cent)^2 + (y_coord(1,ct) - y_cent)^2 );
        end
    end
end
tgt_num = 1 : n_targs;
end
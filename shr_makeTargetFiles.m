function [varargout] = shr_makeTargetFiles(varargin)
% function [varargout] = shr_makeTargetFiles(varargin)
% Creates .tgt files for subject(s) s and block(s) b for the Sequence
% Horizon Reaching experiment
%
% example calls:
%               [G] = shr_makeTargetFiles('s',99, 'b',1, 'use_grid',1, 'grid_size',7);
%               [G] = shr_makeTargetFiles('s',[1:20], 'b',[1:12], 'use_grid',0);
%
% --
% gariani@uwo.ca - 2020.02.14

%% varargin options
s = 99;
b = 1;
use_grid = 1; % choose whether to use a regular grid of predefined targets (1), or not (0)
grid_size = 7; % decide how many crossing lines are splitting the grid, must be odd number
vararginoptions(varargin, {'s', 'b', 'use_grid', 'grid_size'});

%% make target files for all subjects and blocks per subject
G = struct();
for s = s
    S = struct();
    for b = b
        fprintf(1, '\nsubj: %d   block: %d\n', s, b);
        [B] = shr_target(s, b, 'use_grid',use_grid, 'grid_size',grid_size); % B=block (all trials)
        S = addstruct(S, B); % S=subject (all blocks)
    end
    G = addstruct(G, S); % G=group (all subjects)
end
varargout{1} = G;
end

function [varargout] = shr_target(s, b, varargin)
% function [varargout] = shr_target(s, b)
% This function generates .tgt files, one per block
%
% inputs: vector of subject numbers (s), vector of block numbers (b)
% output: saved filename (fn), block structure (B)

%% varargin options
use_grid = 1; % choose whether to use a regular grid of predefined targets (1), or not (0)
grid_size = 7; % decide how many crossing lines are splitting the grid, must be odd number
vararginoptions(varargin, {'use_grid', 'grid_size'});

%% define target folder
tgtDir = '/Users/gariani/Documents/robotcode/projects/SeqHorizonReach/tgt'; %save tgt files in the right folder
if ~exist(tgtDir,'dir'); mkdir(tgtDir); end %create target folder if it doesn't already exist
dtpDir = '/Users/gariani/Documents/robotcode/projects/SeqHorizonReach/dtp'; %save tgt files in the right folder
if ~exist(dtpDir,'dir'); mkdir(dtpDir); end %create target folder if it doesn't already exist

%% experimental details
seq_len = 15; % how many elements per sequence (max)?
n_seq_hor = 4; % how many horizon levels per sequence?

%%% IF YOU CHANGE THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADJUST first entry of <blocktable> in shr_template.dtp %%%%%%%%%%%%%%%
% !! n_trials must be a multiple of n_seq_hor to have a balanced number of
% horizon levels per block (e.g., if n_seq_hor=4, n_trials=4,8,12,16, ...)
n_trials = 40; % how many trials in a block?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_range = [-10 10]; % what's the range of x coords?
y_range = [-10 10]; % what's the range of y coords?
min_d = 22; % what's the minumum distance allowed between successive targets? (expressed in % of largest range)
max_d = 23; % what's the maximum distance allowed between successive targets? (expressed in % of largest range)
if use_grid == 1 % if using the grid, mid and max d are pre-defined according to grid size
    % grid_size MUST BE big enough depending on how many targets in each
    % sequence (seq_len), otherwise there might not be possible solutions
    % to the random search of targets. Also, for how the distance of each
    % reach is currently computed, grid_size MUST BE an odd number.
    grid_size = 7; % decide how many crossing lines are splitting the grid
end

% add sequence horizon factor and randomize order
sh = ceil( (1:n_trials) / (n_trials/n_seq_hor));
sh = sh(randperm(numel(sh)));

%% fill in dataframe structure B for this block
B = struct(); % initialize block structure
trial_mat = []; % initialize empty trial list
for t = 1 : n_trials
    B.TN(t,1) = t; % trial number
    B.BN(t,1) = b; % block number
    B.SN(t,1) = s; % subject number
    % pick a sequence of targets with respective x,y coordinates
    switch use_grid
        case 0
            % either randomly within the predefined range (some constraints apply)
            [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs(seq_len, 'x_range',x_range, 'y_range',y_range, 'min_d',min_d, 'max_d',max_d);
        case 1
            % or randomly but within a predefined grid of potential targets
            [x_coord, y_coord, dist, tgt_num] = pick_nrand_targs_grid(seq_len, 'x_range',x_range, 'y_range',y_range, 'grid_size',grid_size);
    end
    B.seq_len(t,1) = seq_len; % which sequence length?
    B.seq_x(t,:) = x_coord; % what is the x coord of all targets in this sequence?
    B.seq_y(t,:) = y_coord; % what is the x coord of all targets in this sequence?
    B.seq_dist(t,:) = dist; % what is the distance between successive targets in this sequence? (sum is overall path length)
    B.seq_cue(t,:) = tgt_num+1; % which sequence cue? From 2 to n, because 1 is reserved for the home position (see below)
    B.seq_hor(t,1) = sh(t); % which sequence horizon for this sequence?
    B.prep_time_min(t,1)  = 1500; % fixed preparation time (in ms)
    B.prep_time_max(t,1)  = 2000; % fixed preparation time (in ms)
    B.dwell_time_min(t,1) = 150; % fixed dwell time (target acquisition hold) duration (in ms)
    B.dwell_time_max(t,1) = 150; % fixed dwell time (target acquisition hold) duration (in ms)
    B.ITI(t,1)            = 500; % fixed inter-trial-interval duration (in ms)
    B.seq_timeout(t,1)    = 15000; % how much time to complete the whole sequence? (in ms)
    B.tot_trials(t,1)     = n_trials; % how many trials in the entire block?
    B.home(t,1)           = 1; % home target, it's 1 by default (the center of the display)
    this_trial = [B.home(t,1), B.seq_cue(t,:), B.seq_len(t,1), ...
        B.seq_hor(t,1), B.seq_x(t,:), B.seq_y(t,:), B.seq_dist(t,:), ...
        B.prep_time_min(t,1), B.prep_time_max(t,1), ...
        B.dwell_time_min(t,1), B.dwell_time_max(t,1), ...
        B.ITI(t,1), B.seq_timeout(t,1), B.tot_trials(t,1)];
    trial_mat = [trial_mat; this_trial]; %#ok<*AGROW>
    %------------------------------------------------------------------------------------------------------------------------------------
    % sanity check: how does the sequence actually look?
    plot(x_coord',y_coord', 'color',[.88 .88 .88], 'linewidth',3); grid on; hold on;
    plt.scatter(x_coord',y_coord', 'split',tgt_num', 'regression','none', 'label',tgt_num);
    axis equal; xlim(x_range); ylim(y_range); title(sprintf('Horizon: %d', sh(t))); hold off;
    %------------------------------------------------------------------------------------------------------------------------------------
end

%% save structure B as a target file (.tgt) and return output data structure B
outfname = sprintf('shr_s%02d_b%02d', s, b);
dsave(fullfile(tgtDir, sprintf('%s.tgt', outfname)), B);

%% export and save target files in right format for exoskeleton (.dtp)
[dtp] = convert2dtp(trial_mat);
xmlwrite(fullfile(dtpDir, sprintf('%s.dtp', outfname)), dtp);

%% return output data structure B
B.fn = outfname;
varargout{1} = B;
end

function [dtp] = convert2dtp(mat)
%% function [dtp] = convert2dtp(mat)
%
% Convert trial matrix to TP table
mat_size = size(mat);
tp_table = '[';
newline = char(10);
for row=1:mat_size(1)
    tp_table = [tp_table '['];
    for col=1:mat_size(2)
        if col==mat_size(2)
            tp_table = [tp_table num2str(mat(row, col))];
        else
            tp_table = [tp_table num2str(mat(row, col)) ', '];
        end
    end
    if row==mat_size(1)
        tp_table = [tp_table ']'];
    else
        tp_table = [tp_table '],' newline];
    end
end
tp_table = [tp_table ']'];

% read in template .dtp file
dtp = xmlread('rems1_template.dtp');

% replace tp node with new tp_table info
tp_node = dtp.getElementsByTagName('tptable').item(0).item(0);
tp_node.set('Data', tp_table);
end
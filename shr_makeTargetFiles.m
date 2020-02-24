function [varargout] = shr_makeTargetFiles(s, b, varargin)
%
%   This function creates .tgt and .dtp files for the Sequence Horizon
%   Reaching behavioral experiment. In order to function properly, this
%   function requires both the dataframe and the plotlib toolboxes.
%   dataframe:  https://github.com/jdiedrichsen/dataframe.git
%   plotlib:    https://github.com/nejaz1/plotlib.git
%
% INPUTS
%     s:            subject number (can be vector for multiple subjects)
%     b:            block number (can be vector for multiple blocks)
% 
% OPTIONAL VARARGINS
%     n_trials:     defines how many trials per block are created
%     seq_len:      sequence length (number of items / targets)
%     seq_hor:      sequence horizon (how many targets ahead are visible)    
%     x_range:      spatial range of workspace x coordinates [x_min x_max]
%     y_range:      spatial range of workspace y coordinates [y_min y_max]
%     use_grid:     option of whether to use a grid of pre-determined target 
%                   locations (1 = grid locations), or not (0 = random locations)
%     grid_size:    (only matters if using grid) spcifies the resolution of the grid 
%                   (i.e. how many lines are splitting the workspace to create the grid)
%     min_d:        (only matters if NOT using grid) minimum distance between successive 
%                   targets in percentage of largest coord range (default is 20%)
%     max_d:        (only matters if NOT using grid) maximum distance between successive 
%                   targets in percentage of largest coord range (default is 21%)
%     show_seq:     option to plot the sequence or targets for each trial of each block 
%                   (use carefully! e.g., by limiting the number of trials and blocks)
% 
% OUTPUTS
%     G:            output struct for the whole group of subjects [subjects
%                   x blocks x trials] 
%     .tgt files    (one per block per subject) specifying the order and
%                   type of trials (spreadsheet type) 
%     .dtp files    (one per block per subject) specifying information to
%                   be read in by Kinarm Simulink 
%     
% USAGE EXAMPLES
%     [G] = shr_makeTargetFiles(99, 1, 'n_trials',1, 'use_grid',1, 'grid_size',11, 'seq_len',10, 'show_seq',1);
%     [G] = shr_makeTargetFiles(99, 1, 'n_trials',1, 'use_grid',0, 'min_d',5, 'max_d',95, 'seq_len',10, 'show_seq',1);
%     [G] = shr_makeTargetFiles([1:20], [1:12], 'use_grid',0, 'min_d',5, 'max_d',95, 'show_seq',0);
%
% --
% gariani@uwo.ca - 2020.02.24

%% make target files for all subjects and blocks per subject
G = struct();
for s = s
    S = struct();
    for b = b
        fprintf(1, '\nsubj: %d   block: %d\n', s, b);
        [B] = shr_target(s, b, varargin{:}); % B=block (all trials)
        S = addstruct(S, B); % S=subject (all blocks)
    end
    G = addstruct(G, S); % G=group (all subjects)
end
varargout{1} = G;
end

function [varargout] = shr_target(s, b, varargin)
% function [varargout] = shr_target(s, b, varargin)
% This function generates .tgt files, one per block
%
% inputs: vector of subject numbers (s), vector of block numbers (b)
% output: saved filename (fn), block structure (B)

%% define target folder
tgtDir = '/Users/gariani/Documents/robotcode/projects/SeqHorizonReach/tgt'; %save tgt files in the right folder
if ~exist(tgtDir,'dir'); mkdir(tgtDir); end %create target folder if it doesn't already exist
dtpDir = '/Users/gariani/Documents/robotcode/projects/SeqHorizonReach/dtp'; %save tgt files in the right folder
if ~exist(dtpDir,'dir'); mkdir(dtpDir); end %create target folder if it doesn't already exist

%% default experimental details and varargin options
%------------------------------------------------------------------------------------------------------------------------------------
%%% IF YOU CHANGE THIS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ADJUST first entry of <blocktable> in shr_template.dtp %%%%%%%%%%%%%%%
% !! n_trials must be a multiple of seq_hor to have a balanced number of
% horizon levels per block (e.g., if seq_hor=4, n_trials=4,8,12,16, ...)
n_trials = 40; % how many trials in a block?
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------------------------------------------------------------------------------------------------------------------------------
seq_len = 15; % how many elements/targets per sequence (max)?
seq_hor = 4; % how many horizon levels per sequence?
x_range = [-10 10]; % what's the range of x coords?
y_range = [-10 10]; % what's the range of y coords?
min_d = 22; % what's the minumum distance allowed between successive targets? (expressed in % of largest range)
max_d = 23; % what's the maximum distance allowed between successive targets? (expressed in % of largest range)
use_grid = 1; % choose whether to use a regular grid of predefined targets (1), or not (0)
% if using the grid, mid and max d are pre-defined according to grid size
grid_size = 9; % decide how many crossing lines are splitting the grid, must be odd number
% grid_size MUST BE big enough depending on how many targets in each
% sequence (seq_len), otherwise there might not be possible solutions
% to the random search of targets. Also, for how the distance of each
% reach is currently computed, grid_size MUST BE an odd number.
show_seq = 0; % do you want to plot the sequence for each trial in this block? yes (1), or no (0)
vararginoptions(varargin, {'n_trials', 'seq_len', 'seq_hor', 'x_range', 'y_range', 'min_d', 'max_d', 'use_grid', 'grid_size', 'show_seq'});

% add sequence horizon factor and randomize order
sh = ceil( (1:n_trials) / (n_trials/seq_hor));
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
            % or randomly but among a predefined grid of potential targets
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
    if show_seq == 1
        %------------------------------------------------------------------------------------------------------------------------------------
        % sanity check: how does the sequence actually look?
        plot(x_coord',y_coord', 'color',[.88 .88 .88], 'linewidth',3); grid on; hold on;
        plt.scatter(x_coord',y_coord', 'split',tgt_num', 'regression','none', 'label',tgt_num);
        axis image; xlim(x_range); ylim(y_range); title(sprintf('Horizon: %d', sh(t))); hold off;
        set(gcf,'WindowStyle','docked');
        %------------------------------------------------------------------------------------------------------------------------------------
    end
end

%% save structure B as a target file (.tgt) and return output data structure B
outfname = sprintf('shr_s%02d_b%02d', s, b);
dsave(fullfile(tgtDir, sprintf('%s.tgt', outfname)), B);

%% export and save target files in right format for exoskeleton (.dtp)
[dtp] = convert2dtp(trial_mat);
xmlwrite(fullfile(dtpDir, sprintf('%s.dtp', outfname)), dtp);

%% return output structure B
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
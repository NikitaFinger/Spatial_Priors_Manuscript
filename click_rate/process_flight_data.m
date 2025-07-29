%created by Nikita Finger


function click_data = process_flight_data(date_data, trl_num,call_indices_used,trl_name, batID, train_or_test, cond, luminescence, landing, flight_trajectory, release, notes,is_approach_1_present, is_inspect_present, is_approach_2_present, adjust_approach_1_start_frames, adjust_approach_1_end_frames, adjust_inspect_start_frames, adjust_inspect_end_frames, adjust_approach_2_start_frames, adjust_approach_2_end_frames)
addpath 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code'

% Define the full paths to the data files in each folder
bat_pos_file = 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\bat_pos\';
mic_detect_file = 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\mic_detect\';
mic_pos_file = 'Y:\Spatial Memory Project (Xiaoyan and Nikita)\beampattern_analysis\Working Code\mic_pos\';

% Prompt user to select a bat_pos trial file
[filename_bat_pos, path_bat_pos] = uigetfile([bat_pos_file '*' batID '*' date_data '*.mat'], 'Select bat_pos Trial File');
if filename_bat_pos == 0 
    return;
end
load(fullfile(path_bat_pos, filename_bat_pos));

% Prompt user to select a mic_detect trial file
[filename_mic_detect, path_mic_detect] = uigetfile([mic_detect_file '*' batID '*' date_data '*.mat'], 'Select mic_detect Trial File');
if filename_mic_detect == 0
    return;
end
load(fullfile(path_mic_detect, filename_mic_detect));

% Prompt user to select a mic_pos trial file
[filename_mic_pos, path_mic_pos] = uigetfile([mic_pos_file '*' date_data '*.mat'], 'Select mic_pos Trial File');
if filename_mic_pos == 0 
    return;
end
load(fullfile(path_mic_pos, filename_mic_pos));


%% % Create a structure to store the data
click_data = struct();


%% Create a plane 
%creating a plane
frames_w_alldata=find(isfinite(bat_pos{1}(:,1)) & isfinite(bat_pos{2}(:,1))...
  & isfinite(bat_pos{3}(:,1)));

flA=[0 0 0];
flB=[1 0 1];
flC=[-1 0 1];
nnfl=cross(flB-flA,flC-flA);

nn=nan(size(bat_pos{1}));
th_roll=nan(size(bat_pos{1},1),1);
head_vec=nan(size(bat_pos{1}));
for fr=frames_w_alldata'
  %three points on the triangle
  A=bat_pos{1}(fr,:);
  B=bat_pos{2}(fr,:);
  C=bat_pos{3}(fr,:);
  
  %get roll vector from plane
  nn(fr,:)=cross(B-A,C-A);
  th_roll(fr)=atan2(norm(cross(nn(fr,:),nnfl)),dot(nn(fr,:),nnfl));
  
  %direction vector
  midptBC=(B+C)./2;
  head_vec(fr,:)=A-midptBC;
end
%% 
% Extract the xyz coordinates from the cell matrix based on availability of data
if any(isfinite(bat_pos{1}(:, 1)))
    xyz = bat_pos{1}(:, 1:3);
elseif any(isfinite(bat_pos{2}(:, 1)))
    xyz = bat_pos{2}(:, 1:3);
elseif any(isfinite(bat_pos{3}(:, 1)))
    xyz = bat_pos{3}(:, 1:3);
else
    error('No valid data found in bat_pos cell matrix.');
end


% Extrapolate missing XYZ coordinates using linear interpolation
for i = 1:3 % Loop over X, Y, and Z dimensions
    nan_indices = find(isnan(xyz(:, i))); % Find NaN indices in the current dimension
    if ~isempty(nan_indices)
        for j = 1:length(nan_indices)
            nan_idx = nan_indices(j);
            before_nan_idx = find(isfinite(xyz(1:nan_idx, i)), 1, 'last');
            after_nan_idx = find(isfinite(xyz(nan_idx:end, i)), 1);
            if ~isempty(before_nan_idx) && ~isempty(after_nan_idx)
                % Linear interpolation
                interpolated_value = interp1([before_nan_idx, after_nan_idx + nan_idx - 1], ...
                    [xyz(before_nan_idx, i), xyz(after_nan_idx + nan_idx - 1, i)], nan_idx);
                xyz(nan_idx, i) = interpolated_value;
            end
        end
    end
end
call_locs = cat(1, call.locs);
call_locs = call_locs(:, 1) / 2500;

% Define the position of the perch
perch_x = mic_pos(51, 1); % Make sure this is the correct index for the perch position

%% 

% Calculate flight speed
speed_threshold = 0.1; % Set your desired speed threshold
flight_speed = sqrt(diff(xyz(:, 1)).^2 + diff(xyz(:, 2)).^2 + diff(xyz(:, 3)).^2);

% Find the index where flight starts based on speed
start_flight_index = find(flight_speed > speed_threshold, 1, 'first');

% Find the index where flight ends based on speed
end_flight_index = find(flight_speed > speed_threshold, 1, 'last');

% Define the start position based on flight start index
start_x = xyz(start_flight_index, 1);
start_y = xyz(start_flight_index, 2);
start_z = xyz(start_flight_index, 3);

% Define the end position based on flight end index
end_position_x = xyz(end_flight_index, 1);
end_position_y = xyz(end_flight_index, 2);
end_position_z = xyz(end_flight_index, 3);

% Create a new figure in the live script output
f = figure('Name', 'Flight Trajectory', 'NumberTitle', 'off');

% Plot the trajectory of the bat starting from the green start point
plot3(xyz(:, 1), xyz(:, 2), xyz(:, 3), '-b'); % Trajectory of the bat in blue
hold on
scatter3(perch_x, mic_pos(51, 2), mic_pos(51, 3), 100, 'r', 'filled', 'MarkerEdgeColor', 'k');
grid on;
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Flight Trajectory');
legend('Trajectory', 'Perch');
grid on
axis equal
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');

% Mark the start and end of the trajectory
scatter3(start_x, start_y, start_z, 100, 'g', 'filled', 'MarkerEdgeColor', 'k');
scatter3(end_position_x, end_position_y, end_position_z, 100, 'b', 'filled', 'MarkerEdgeColor', 'k');
text(start_x, start_y, start_z, 'Start', 'Color', 'g', 'FontSize', 12, 'FontWeight', 'bold');
text(end_position_x, end_position_y, end_position_z, 'End', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');

% Display the figure in the live script output
display(f);

% Calculate the distances between the bat positions and the perch using the
% Euclideans distance forumula considering all three dimensions (x,y,z). 
distances_to_perch = sqrt(sum((xyz - repmat([perch_x, mic_pos(51, 2), mic_pos(51, 3)], size(xyz, 1), 1)) .^ 2, 2));

%% Approach 1

if is_approach_1_present
% Find the Approach_1_start using the time when the bat begins to fly (first flight frame with xyz coordinate)
approach_1_start = find(isfinite(xyz(:, 1)), 1) + adjust_approach_1_start_frames; 
% Define the time threshold for distance increase detection
time_threshold = 0.03; % Default time threshold of 10 milliseconds (assuming each frame corresponds to 1 millisecond)

% Find the Approach_1_end by finding the index where the distance starts to increase again for at least the specified time threshold
increasing_distance_mask = [false; diff(distances_to_perch) > 0];
approach_1_end_indices = find(increasing_distance_mask(approach_1_start:end) & (cumsum(increasing_distance_mask(approach_1_start:end)) * time_threshold >= 1.5), 1);
approach_1_end_indices = approach_1_end_indices - adjust_approach_1_end_frames 
if ~isempty(approach_1_end_indices)
    % If increase found, set the approach_1_end as the first frame of increase
    approach_1_end = approach_1_start + approach_1_end_indices - 1;
else
    % If no increase found, find the first position the bat lands on or gets to the perch
    approach_1_end = find(distances_to_perch(approach_1_start:end) == min(distances_to_perch(approach_1_start:end)), 1);
    
    if ~isempty(approach_1_end)
        approach_1_end = approach_1_start + approach_1_end - 1;
    else
        % If no such position found, set end as the last position
        approach_1_end = size(xyz, 1);
    end
end

% Find the time duration for the first approach
% Switched approach_1_start with approach_1_end
duration_approach_1_s = approach_1_end - approach_1_start + 1; % Duration of approach to the perch 1 in seconds
duration_approach_1_s = duration_approach_1_s/100; 

% Find the indices of calls that fall within the time duration of the first approach
calls_approach_1_indices = find(call_locs >= approach_1_start & call_locs <= approach_1_end);

% Calculate the click rate for Approach 1
click_rate_approach_1 = numel(calls_approach_1_indices) / duration_approach_1_s; % Click rate in Approach 1 in clicks per second

% Mark the Approach_1 section on the plot in black color
plot3(xyz(approach_1_start:approach_1_end, 1), xyz(approach_1_start:approach_1_end, 2), xyz(approach_1_start:approach_1_end, 3), 'k', 'LineWidth', 2, 'DisplayName', 'Start');

% Save Approach_1 start and end points for later use if needed
approach_1_start_x = xyz(approach_1_start, 1);
approach_1_start_y = xyz(approach_1_start, 2);
approach_1_start_z = xyz(approach_1_start, 3);

approach_1_end_x = xyz(approach_1_end, 1);
approach_1_end_y = xyz(approach_1_end, 2);
approach_1_end_z = xyz(approach_1_end, 3);

% Define the XYZ coordinate for the end of Approach 1 phase
xyz_approach_1_end = xyz(approach_1_end, :);


% Label the start points of Approach_1
text(approach_1_start_x, approach_1_start_y, approach_1_start_z, 'start', 'Color', 'b', 'FontSize', 12, 'FontWeight', 'bold');


% Add the legend entry for Approach_1
legend('Trajectory', 'Perch', 'Approach 1');
fprintf('Click Rate during Approach 1: %.4f clicks per second\n', click_rate_approach_1);

% Calculate the time from each click to the end of Approach 1
clicks_approach_1_times_to_end = (approach_1_end - call_locs(calls_approach_1_indices)) / frame_rate;

% Calculate the distance to the perch for each click in Approach 1
clicks_approach_1_distances_to_perch = distances_to_perch(calls_approach_1_indices);

% Store these values in the click_data structure for Approach 1
click_data.approach_1_click_times_to_end = clicks_approach_1_times_to_end;
click_data.approach_1_click_distances_to_perch = clicks_approach_1_distances_to_perch;


% % Calculate distances from the end of Approach 1 phase when a call is made
approach_1_end_xyz = xyz(approach_1_end, :);
call_locs_rounded = round(call_locs(calls_approach_1_indices));
distances_from_end_approach_1 = sqrt(sum((approach_1_end_xyz - xyz(call_locs_rounded, :)).^2, 2));
%calculate number of call in approach 1 phase 
num_calls_approach_1 = numel(calls_approach_1_indices);

%Store click rate data
click_data.click_rate_approach_1 = click_rate_approach_1;

% Store phase durations
click_data.duration_approach_1_s = duration_approach_1_s;

%Store number of calls
click_data.num_calls_approach_1 = num_calls_approach_1;

% Store call indices and timestamps for each phase
click_data.calls_approach_1_indices = calls_approach_1_indices;
click_data.call_locs_approach_1 = call_locs(calls_approach_1_indices);

%  Store distances_from_end_approach_1  in the click_data structure
 click_data.distances_from_end_approach_1 = distances_from_end_approach_1;
end 
%% Inspect 

 if is_inspect_present   
% Find the Inspect_start (when the bat passes the perch on the Y-plane the first time)
inspect_start = find(xyz(approach_1_end:end, 2) > mic_pos(51, 2), 1) + adjust_inspect_start_frames + approach_1_end - 1;


% Find the Inspect_end (when the bat passes the perch on the Y-plane the second time)
% Find the end of the Inspect phase
inspect_end = find(xyz(inspect_start:end, 2) < mic_pos(51, 2), 1) + adjust_inspect_end_frames + inspect_start - 1;


% Find the time duration for the inspect phase
% Switched inspect_start with inspect_end
duration_inspect_s = inspect_end - inspect_start + 1; % Duration of approach to the perch 1 in seconds
duration_inspect_s = duration_inspect_s/100; 

% Find the indices of calls that fall within the time duration of the inspect phase
calls_inspect_indices = find(call_locs >= inspect_start & call_locs <= inspect_end);

% Calculate the click rate for the explore phase
click_rate_inspect = numel(calls_inspect_indices) / duration_inspect_s; % Click rate in Explore phase in clicks per second

% Mark the Inspect section on the plot in purple color
plot3(xyz(inspect_start:inspect_end, 1), xyz(inspect_start:inspect_end, 2), xyz(inspect_start:inspect_end, 3), '-m', 'LineWidth', 2);

fprintf('Click Rate during Inspect: %.4f clicks per second\n', click_rate_inspect);
% Add the legend entry for Inspect
legend('Trajectory', 'Perch', 'Approach 1', 'Inspect');
 
 % Calculate time from the end of each phase
 time_from_end_inspect = min(duration_inspect_s, (call_locs(calls_inspect_indices) - inspect_end) / 100);
 
 %Calculate number of calls in inspect phase
num_calls_inspect = numel(calls_inspect_indices);

%Store click rate data
click_data.click_rate_inspect = click_rate_inspect;

% Store phase durations
click_data.duration_inspect_s = duration_inspect_s;

%Store number of calls
click_data.num_calls_inspect = num_calls_inspect;

% Store call indices and timestamps for each phase
click_data.calls_inspect_indices = calls_inspect_indices;
click_data.call_locs_inspect = call_locs(calls_inspect_indices);

 end



%% Approach 2 

if is_approach_2_present
    
if is_inspect_present == 1
    % Find Approach_2_start (when the bat starts flying towards the perch again after the inspect phase)
    approach_2_start_candidates = find(diff(distances_to_perch(inspect_end:end)) < 0); % Find instances where the distance to perch decreases

    if isempty(approach_2_start_candidates)
        approach_2_start = min(inspect_end, end_position);
    else
        approach_2_start = inspect_end + approach_2_start_candidates(1) - 1 + adjust_approach_2_start_frames; % Use the first instance of decreasing distance as the start of Approach 2
    end
else
    % If there is no inspect phase, set Approach 2 start based on distance decrease
    approach_2_start_candidates = find(diff(distances_to_perch(approach_1_end:end)) < 0); % Find instances where the distance to perch decreases

    if isempty(approach_2_start_candidates)
        approach_2_start = min(approach_1_end, end_position);
    else
        approach_2_start = approach_1_end + approach_2_start_candidates(1) - 1 + adjust_approach_2_start_frames; % Use the first instance of decreasing distance as the start of Approach 2
    end
end

% Find Approach_2_end (when the bat passes the perch on the Y-plane the third time)
approach_2_end_candidates = find(xyz(approach_2_start:end, 2) < mic_pos(51, 2)); % Find all instances the bat goes below the perch

if isempty(approach_2_end_candidates)
    approach_2_end = min(approach_1_end, end_position);
else
    approach_2_end = approach_2_start + approach_2_end_candidates(end) - 1 + adjust_approach_2_end_frames; % Use the last instance of going below the perch as the end of Approach 2


% Find the time duration for the second approach
duration_approach_2_s = approach_2_end - approach_2_start + 1; % Duration of approach 2
duration_approach_2_s = duration_approach_2_s/100; 

% Find the indices of calls that fall within the time duration of approach 2 phase
calls_approach_2_indices = find(call_locs >= approach_2_start & call_locs <= approach_2_end);

% Calculate the click rate for the approach 2 phase
click_rate_approach_2 = numel(calls_approach_2_indices) / duration_approach_2_s; % Click rate in Explore phase in clicks per second


% Plot Approach 2 section on the plot in magenta color
plot3(xyz(approach_2_start:approach_2_end, 1), xyz(approach_2_start:approach_2_end, 2), xyz(approach_2_start:approach_2_end, 3), 'g', 'LineWidth', 2);

fprintf('Click Rate during Approach 2: %.4f clicks per second\n', click_rate_approach_2);

% Calculate time from the end of each phase
time_from_end_approach_2 = min(duration_approach_2_s, (call_locs(calls_approach_2_indices) - approach_2_end) / 100);

% % Calculate distances from the end of Approach 2 phase when a call is made
approach_2_end_xyz = xyz(approach_2_end, :);
call_locs_rounded = round(call_locs(calls_approach_2_indices));
distances_from_end_approach_2 = sqrt(sum((approach_2_end_xyz - xyz(call_locs_rounded, :)).^2, 2));

 %Calculate number of calls in approach 2 phase 
 num_calls_approach_2 = numel(calls_approach_2_indices);
 
 %Store approach 2 click rate
 click_data.click_rate_approach_2 = click_rate_approach_2;
 
 % Store phase durations
 click_data.duration_approach_2_s = duration_approach_2_s;
 
 %Store number of calls
click_data.num_calls_approach_2 = num_calls_approach_2;

% Store call indices and timestamps for each phase
click_data.calls_approach_2_indices = calls_approach_2_indices;
click_data.call_locs_approach_2 = call_locs(calls_approach_2_indices);

%  Store distances_from_end_approach_2  in the click_data structure
click_data.distances_from_end_approach_2 = distances_from_end_approach_2;
    end
end

% Only proceed if the Approach 2 phase is present
if is_approach_2_present
    % Calculate the time from each click to the end of Approach 2
    clicks_approach_2_times_to_end = (approach_2_end - call_locs(calls_approach_2_indices)) / frame_rate;
    
    % Calculate the distance to the perch for each click in Approach 2
    clicks_approach_2_distances_to_perch = distances_to_perch(calls_approach_2_indices);
    
    % Store these values in the click_data structure for Approach 2
    click_data.approach_2_click_times_to_end = clicks_approach_2_times_to_end;
    click_data.approach_2_click_distances_to_perch = clicks_approach_2_distances_to_perch;
end
% % %% 
%% Store Index Values
click_data.batID = batID; %as string
click_data.cond = cond; % as string
click_data.date_data = date_data ; %as string
click_data.trl_num = trl_num'; %as string, vicon 
click_data.train_or_test = train_or_test;  %as string
click_data.cond = cond; %as string, experiment condition
click_data.luminescence = luminescence; %as string
click_data.landing = landing; %as string
click_data.flight_trajectory =flight_trajectory; %as string
click_data.release = release; %as string
click_data.notes =notes;% as string
click_data.trl_name = trl_name; % as string


click_data.is_approach_1_present = is_approach_1_present;  % as string
click_data.is_inspect_present = is_inspect_present;     % as string
click_data.is_approach_2_present = is_approach_2_present;% as string


%% 

% Plot the clicks as blue sticks aligned with the bat's trajectory
stick_length = 0.1; % Length of the sticks representing the calls    
   


% Check if Inspect and Approach 2 phase data are present
is_inspect_present = isfield(click_data, 'click_rate_inspect');
is_approach_2_present = isfield(click_data, 'click_rate_approach_2');

% Create a legend_labels cell array with the default entries
legend_entries = {'Trajectory', 'Perch', 'Approach 1'};

% Add "Inspect" to the legend_labels if it's present
if is_inspect_present
    legend_entries = [legend_entries, 'Inspect'];
end

% Add "Approach 2" to the legend_labels if it's present
if is_approach_2_present
    legend_entries = [legend_entries, 'Approach 2'];
end

% Add "Click" to the legend_labels
% legend_entries = [legend_entries, 'Click'];

% Plot the legend using the legend_labels
legend(legend_entries, 'Location', 'Best');

hold on;

% Plot the calls as sticks
for i = 1:length(call_locs)
    frame_number = round(call_locs(i));
    if frame_number >= 1 && frame_number <= length(xyz) % Use length(xyz) instead of length(x)
        call_x = xyz(frame_number, 1); % Use xyz(:, 1) for X-coordinates
        call_y = xyz(frame_number, 2); % Use xyz(:, 2) for Y-coordinates
        call_z = xyz(frame_number, 3); % Use xyz(:, 3) for Z-coordinates
        if i == 1
            % Plot the first call with legend entry
            plot3([call_x, call_x], [call_y, call_y], [call_z, call_z + stick_length], '-r', 'LineWidth', 2, 'DisplayName', 'Clicks');
        else
            % Plot subsequent calls without legend entry
            plot3([call_x, call_x], [call_y, call_y], [call_z, call_z + stick_length], '-r', 'LineWidth', 2);
        end
    end
end

%% Extract entire distances of each approach 
% Calculate the entire distance of Approach 1
distance_approach_1 = sum(sqrt(sum(diff(xyz(approach_1_start:approach_1_end, :)).^2, 2)));

% Calculate the entire distance of Inspect phase
if is_inspect_present
    distance_inspect = sum(sqrt(sum(diff(xyz(inspect_start:inspect_end, :)).^2, 2)));
else
    distance_inspect = [];
end

% Calculate the entire distance of Approach 2
if is_approach_2_present
    distance_approach_2 = sum(sqrt(sum(diff(xyz(approach_2_start:approach_2_end, :)).^2, 2)));
else
    distance_approach_2 = [];
end




% Calculate total time traveled
total_time_traveled_s = (end_flight_index - start_flight_index) / 1000; % Total time traveled in seconds

% Calculate total distance traveled
% Calculate total distance traveled
total_distance_traveled_m = sqrt((end_position_x - start_x) + (end_position_y - start_y) + (end_position_z - start_z)); % Total distance traveled in meters
 % Total distance traveled in meters; % Total distance traveled in meters


% Store distances and times in the click_data structure
click_data.distance_approach_1 = distance_approach_1;
click_data.distance_inspect = distance_inspect;
click_data.distance_approach_2 = distance_approach_2;
click_data.total_distance_traveled_m = total_distance_traveled_m;

click_data.total_time_traveled_s = total_time_traveled_s;


%% % Calculate click rates
total_duration_s = 8; % Total duration in seconds
click_rate_total = numel(call_locs) / total_duration_s; % Click rate in clicks per second

% Assuming a constant frame rate, e.g., 250 frames per second
frame_rate = 250; 

% Convert call locations to time
call_times = call_locs / frame_rate;

% Check if Approach 1 is present and if calls are within Approach 1
  % Find calls within Approach 1
    calls_in_approach_1 = call_locs >= approach_1_start & call_locs <= approach_1_end;
    
    % Extract the call locations within Approach 1
    call_locs_approach_1 = call_locs(calls_in_approach_1);
    
    % Calculate time from each click to the end of Approach 1
    time_to_end_approach_1 = (approach_1_end - call_locs_approach_1) / frame_rate;

if is_approach_1_present
    % Logical array for calls within Approach 1
    calls_in_approach_1 = call_locs >= approach_1_start & call_locs <= approach_1_end;
    
    % Extract call locations within Approach 1
    call_locs_approach_1 = call_locs(calls_in_approach_1);
    
    % Calculate time from each click to the end of Approach 1
    time_to_end_approach_1 = (approach_1_end - call_locs_approach_1) / frame_rate;
    
    % Find the indices of the calls within Approach 1. Since calls_in_approach_1 is a logical array,
    % use find() to convert it to indices.
    indices_approach_1 = find(calls_in_approach_1);
  
    click_data.approach_1_click_indices_and_times = [indices_approach_1, time_to_end_approach_1];
end



% Check if Approach 2 is present and if calls are within Approach 2
if is_approach_2_present
    % Logical array for calls within Approach 2
    calls_in_approach_2 = call_locs >= approach_2_start & call_locs <= approach_2_end;
    
    % Extract call locations within Approach 2
    call_locs_approach_2 = call_locs(calls_in_approach_2);
    
    % Calculate time from each click to the end of Approach 2
    time_to_end_approach_2 = (approach_2_end - call_locs_approach_2) / frame_rate;
    
    % Convert the logical array to indices
    indices_approach_2 = find(calls_in_approach_2);
    
    % Store times to end and the corresponding indices for Approach 2 calls in the click_data structure.
    click_data.approach_2_click_indices_and_times = [indices_approach_2, time_to_end_approach_2];
end


% Initialize call pairs array
call_pairs = [];


%% Save data 

% Define the base directory
base_directory = 'Y:\\Spatial Memory Project (Xiaoyan and Nikita)\\Analysis\\Click Rate\\Analysis\\mat files';

% Create the file name using trl_name
file_name = [trl_name, '_click_data.mat'];  % Assuming trl_name is a string

% Combine the directory and file name to create the full file path
full_save_path = fullfile(base_directory, file_name);

% Save the click_data structure to the specified path, overwriting any existing file
save(full_save_path, 'click_data', '-v7.3')

fprintf('Total Click Rate: %.4f clicks per second\n', click_rate_total);
% Display the click rates


end






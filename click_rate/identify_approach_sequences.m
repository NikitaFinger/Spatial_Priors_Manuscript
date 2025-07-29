function identify_approach_sequences(date_data, trl_num, batID, train_or_test, cond, luminescence, landing, flight_trajectory, release, notes, is_approach_1_present, is_inspect_present, is_approach_2_present, adjust_approach_1_start_frames, adjust_approach_1_end_frames, adjust_inspect_start_frames, adjust_inspect_end_frames, adjust_approach_2_start_frames, adjust_approach_2_end_frames)
%% Approach 1

if is_approach_1_present
% Find the Approach_1_start using the time when the bat begins to fly (first flight frame with xyz coordinate)
approach_1_start = find(isfinite(xyz(:, 1)), 1) + adjust_approach_1_start_frames; 
% Define the time threshold for distance increase detection
time_threshold = 0.01; % Default time threshold of 10 milliseconds (assuming each frame corresponds to 1 millisecond)

% Find the Approach_1_end by finding the index where the distance starts to increase again for at least the specified time threshold
increasing_distance_mask = [false; diff(distances_to_perch) > 0];
approach_1_end_indices = find(increasing_distance_mask(approach_1_start:end) & (cumsum(increasing_distance_mask(approach_1_start:end)) * time_threshold >= 0.0), 1);

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


%% Inspect 

 if is_inspect_present   
% Find the Inspect_start (when the bat passes the perch on the Y-plane the first time)
inspect_start = find(xyz(approach_1_end:end, 2) > mic_pos(51, 2), 1) + adjust_inspect_start_frames + approach_1_end - 1;


% Find the Inspect_end (when the bat passes the perch on the Y-plane the second time)
% Find the end of the Inspect phase
inspect_end = find(xyz(inspect_start:end, 2) < mic_pos(51, 2), 1) + adjust_inspect_end_frames + inspect_start - 1;




% Find the indices of calls that fall within the time duration of the inspect phase
calls_inspect_indices = find(call_locs >= inspect_start & call_locs <= inspect_end);


% Mark the Inspect section on the plot in purple color
plot3(xyz(inspect_start:inspect_end, 1), xyz(inspect_start:inspect_end, 2), xyz(inspect_start:inspect_end, 3), '-m', 'LineWidth', 2);





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



% Find the indices of calls that fall within the time duration of approach 2 phase
calls_approach_2_indices = find(call_locs >= approach_2_start & call_locs <= approach_2_end);



% Plot Approach 2 section on the plot in magenta color
plot3(xyz(approach_2_start:approach_2_end, 1), xyz(approach_2_start:approach_2_end, 2), xyz(approach_2_start:approach_2_end, 3), 'g', 'LineWidth', 2);


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
legend_entries = [legend_entries, 'Click'];

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


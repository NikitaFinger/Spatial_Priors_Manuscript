%% % Calculate click rates
total_duration_s = 8; % Total duration in seconds
click_rate_total = numel(call_locs) / total_duration_s; % Click rate in clicks per second

%% % Create a structure to store the data
click_data = struct();


% Display the click rates
fprintf('Total Click Rate: %.4f clicks per second\n', click_rate_total);%Approach 1 
% Calculate the click rate for Approach 1
click_rate_approach_1 = numel(calls_approach_1_indices) / duration_approach_1_s; % Click rate in Approach 1 in clicks per second

% Calculate time from the end of each phase
time_from_end_approach_1 = min(duration_approach_1_s, (call_locs(calls_approach_1_indices) - approach_1_end) / 100);

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
 
 % Find the time duration for the inspect phase
% Switched inspect_start with inspect_end
duration_inspect_s = inspect_end - inspect_start + 1; % Duration of approach to the perch 1 in seconds
duration_inspect_s = duration_inspect_s/100; 

% Calculate the click rate for the explore phase
click_rate_inspect = numel(calls_inspect_indices) / duration_inspect_s; % Click rate in Explore phase in clicks per second

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
 
%% %approach 2
% Find the time duration for the second approach
duration_approach_2_s = approach_2_end - approach_2_start + 1; % Duration of approach 2
duration_approach_2_s = duration_approach_2_s/100; 

% Calculate the click rate for the approach 2 phase
click_rate_approach_2 = numel(calls_approach_2_indices) / duration_approach_2_s; % Click rate in Explore phase in clicks per second

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

 
%This code was developed to open pre-stored OCT phase data from a needle 
%probe and calculate the phase difference between subsequent A-lines. 
%Smoothing is applied. 

%V1.0, Danielle J. Harper, Yongjoo Kim, Alejandra Gómez-Ramírez and Benjamin J. Vakoc, 2023
%Further information can be found in this paper:
   % https://doi.org/10.48550/arXiv.2305.14390
% Set this to 1 for example data set or 2 for salmon data
% Clear the environment and set up for opening the phase files
clc;
close all;
clear all;

case_file = 'UW_pork_belly'; 

% Switch between cases
switch case_file
    case 'example_salmon' % Example data set       
        % Set directory and filenames for phase files
        data_files_folder = 'E:\My Drive\4_UW-Madison\2_4_MGH+FDA+UW OCT\needle_code_share\Example_data';
        file_name1 = '[p.needle_2][s.salmon][06-16-2021_14-15-13]phase1.mgh';
        file_name2 = '[p.needle_2][s.salmon][06-16-2021_14-15-13]phase2.mgh';
        read_opt.nFrames = 1024;
        machine_ID = 'MGH';
        % Parameters for Region of Interest
        roi_start_row = 220; % Starting row of the region of interest
        roi_end_row = 272;   % Ending row of the region of interest
        
        % Optional: Parameters for Calibration Line (if needed)
        calibration_line_start = 205; % Start of calibration line (fiber tip)
        calibration_line_end = 215;   % End of calibration line (fiber tip)

               
    case 'UW_salmon' % Salmon data
        % Set directory and filenames for phase files
        data_files_folder = 'H:\OFDIData\user.Ricardo\[p.231220_PS_Needle_Probe]\[p.231220_PS_Needle_Probe][s.Salmon_Probe_04_Test3_Dist1.2cm][12-20-2023_15-52-35]';
        file_name1 = '[p.231220_PS_Needle_Probe][s.Salmon_Probe_04_Test3_Dist1.2cm][12-20-2023_15-52-35].phaseXA.mgh';
        file_name2 = '[p.231220_PS_Needle_Probe][s.Salmon_Probe_04_Test3_Dist1.2cm][12-20-2023_15-52-35].phaseXB.mgh';
        read_opt.nFrames = 1023;
        machine_ID = 'SPARC';
        % Parameters for Region of Interest
        roi_start_row = 476; % Starting row of the region of interest
        roi_end_row = 526;   % Ending row of the region of interest
        
        % Optional: Parameters for Calibration Line (if needed)
        calibration_line_start = 466; % Start of calibration line (fiber tip)
        calibration_line_end = 476;   % End of calibration line (fiber tip)

    case 'UW_pork_belly' %Pork Belly data collected 231227
        % Set directory and filenames for phase files
        data_files_folder = 'H:\OFDIData\user.Ricardo\[p.240108_Needle_Probe_Data_Set1]\[p.231227_PS_Needle_Porcine][s.Belly_Probe_5_Test_1][12-27-2023_14-06-57]';
        file_name1 = '[p.231227_PS_Needle_Porcine][s.Belly_Probe_5_Test_1][12-27-2023_14-06-57].phaseXA.mgh';
        file_name2 = '[p.231227_PS_Needle_Porcine][s.Belly_Probe_5_Test_1][12-27-2023_14-06-57].phaseXB.mgh';
        read_opt.nFrames = 1023;
        machine_ID = 'SPARC';
        % Parameters for Region of Interest
        roi_start_row = 548; % Starting row of the region of interest
        roi_end_row = 598;   % Ending row of the region of interest
        
        % Optional: Parameters for Calibration Line (if needed)
        calibration_line_start = 538; % Start of calibration line (fiber tip)
        calibration_line_end = 548;   % End of calibration line (fiber tip)
        
    otherwise
        error('Invalid case name. Please specify `example` for example data set or `salmon` for salmon data.');
end
addpath(data_files_folder);

% Setting the directory for reading the file
read_opt.dirname = data_files_folder;

% Define initial frame and number of frames to read
read_opt.iFrame = 1;

% Read the MGH file with the specified options
phase_img_8bit1 = readMgh(file_name1, read_opt);

%% Convert phase data to true phase values and calculate phase difference

% Convert the 8-bit phase images to true phase values in the range of -pi to pi.
% The original data was saved as 8-bit (0-255) and is now being converted back.
number_bits = 8;
phase1 = single(phase_img_8bit1).*(2*pi)./(2^number_bits-1)-pi;

% Retrieve the dimensions of the phase image (rows, columns, frames)
[~, ~, num_frm] = size(phase1); 

% Define column distance based on machine type
col_distance = (strcmp(machine_ID, 'MGH')) * 2 + (strcmp(machine_ID, 'SPARC')) * 2;

% Call the function to calculate phase difference
diff_ph1 = calculate_phase_difference(phase1, col_distance);


%% Slice the phase difference array to focus on a specific region of interest by 
% phase difference of only the region of interest (sliced)

% Calculate the number of rows and columns based on the region of interest
num_rows = roi_end_row - roi_start_row + 1; % Number of rows
num_rows_cal = calibration_line_end - calibration_line_start; % Number of rows for calibration

% Slice the phase difference array to focus on the specified region of interest
diff_ph1_sliced = diff_ph1(roi_start_row:roi_end_row,:,:); % Slicing


%% Concatenate two frames together, handling cases with an odd number of frames
diff_ph_con = concatenate_frames(diff_ph1_sliced);
checksum = sum(diff_ph_con(:));
disp(['Checksum after concatenating: ', num2str(checksum)]);
% Calculating the number of frames after concatenation
num_frm_con = ceil(num_frm/2);
%% Apply thresoldingg
% Parameters for Filtering
threshold_value = 2.5; % Threshold value for filtering
filtered_diff_ph_con = diff_ph_con;
filtered_diff_ph_con(filtered_diff_ph_con > threshold_value) = 0;
num_modifications = sum(diff_ph_con(:) > threshold_value);

% Display the number of modifications (if needed)
disp(['Number of modifications made during filtering: ', num2str(num_modifications)]);

checksum = sum(filtered_diff_ph_con(:));
disp(['Checksum after thresolding: ', num2str(checksum)]);
%% apply Moving median filters 

% Parameters
medianFilterSize = 73; % Size of the window for the moving median e.g. 73
edgeExtensionSize = 36; % Edge extension size for concatenation e.g. 36
columnSize = size(filtered_diff_ph_con,2); % Defined elsewhere in your script

% Initialize arrays with the parameterized sizes
movm = zeros(num_rows, columnSize);
movma = zeros(num_rows, columnSize + edgeExtensionSize);
movmb = zeros(num_rows, columnSize + 2 * edgeExtensionSize);
diff_phase_con_med = zeros(num_rows, columnSize, num_frm_con);
diff_phase_con_med_a = zeros(num_rows, columnSize + edgeExtensionSize, num_frm_con);
diff_phase_con_med_b = zeros(num_rows, columnSize + 2 * edgeExtensionSize, num_frm_con);

% Processing loop
for k = 1:num_frm_con
    if num_frm_con == 1
        movm = movmedian(filtered_diff_ph_con(:,:,k), medianFilterSize, 2);
        diff_phase_con_med(:,:,k) = movm;
        break;
    end
    if k + 1 > num_frm_con
        diff_phase_con_med_a(:,:,k) = horzcat(filtered_diff_ph_con(:,columnSize-edgeExtensionSize+1:columnSize,k-1), filtered_diff_ph_con(:,:,k));
        movma = movmedian(diff_phase_con_med_a(:,:,k), medianFilterSize, 2);
        diff_phase_con_med(:,:,k) = movma(:,edgeExtensionSize+1:end);
        break;
    end

    if k ~= 1
        diff_phase_con_med_b(:,:,k) = horzcat(filtered_diff_ph_con(:,columnSize-edgeExtensionSize+1:columnSize,k-1), filtered_diff_ph_con(:,:,k), filtered_diff_ph_con(:,1:edgeExtensionSize,k+1));
        movmb = movmedian(diff_phase_con_med_b(:,:,k), medianFilterSize, 2);
        diff_phase_con_med(:,:,k) = movmb(:,edgeExtensionSize+1:columnSize+edgeExtensionSize);
        continue;
    end

    diff_phase_con_med_a(:,:,k) = horzcat(filtered_diff_ph_con(:,:,k), filtered_diff_ph_con(:,1:edgeExtensionSize,k+1));
    movma = movmedian(diff_phase_con_med_a(:,:,k), medianFilterSize, 2);
    diff_phase_con_med(:,:,k) = movma(:,1:columnSize,:);

    

end

checksum = sum(diff_phase_con_med(:));
disp(['Checksum after movmedian: ', num2str(checksum)]);

% Add pi to the filtered phase difference for display
diff_phase_con_med_shifted = diff_phase_con_med + pi;
figure('Name', 'movmedian+2pi'), imshow3D(diff_phase_con_med_shifted, [0 2*pi]);
colormap(redblue);
%% Average across B scans
% Compute the average of each A-scan in the region of interest for each frame
num_cols = size(filtered_diff_ph_con, 2); % Number of columns in the filtered_diff_ph_con array
avg_a_scan_per_frame = zeros(1, num_cols, num_frm_con);

for frame_idx = 1:num_frm_con
    % Compute the mean along the first dimension (rows)
    avg_a_scan_per_frame(:,:,frame_idx) = mean(diff_phase_con_med(:,:,frame_idx), 1);
end

% Concatenate the averages from all frames into a single row
concatenated_avg_a_scan = reshape(avg_a_scan_per_frame, [1, num_cols * num_frm_con]);

% Normalize by subtracting the mean across the concatenated data
normalized_avg_a_scan = concatenated_avg_a_scan - mean(concatenated_avg_a_scan, 2);
checksum = sum(normalized_avg_a_scan(:));
disp(['Checksum after Averaging across B scans and more: ', num2str(checksum)]);

% Plotting the unwrapped concatenated and normalized averages
figure('Name', 'Average A-Scan per Frame');
plot(unwrap(normalized_avg_a_scan(:,:)));
title('Average A-Scan Across Frames');
xlabel('Concatenated Frame Index');
ylabel('Normalized Average Value');
%% Calculate Distance

% Parameters
wavelength = 0.0013; % Wavelength in mm (1.3 um)
refractive_index = 1.39; % Refractive index of tissue

% Applying the formula to convert normalized phase change to distance
% based on wavelength and refractive index
distance = (normalized_avg_a_scan * wavelength) / (4 * pi * refractive_index);

checksum = sum(distance(:));
disp(['Checksum after Calculating distance: ', num2str(checksum)]);

% Correct for any drift if necessary
drift_correction = 1.6e-7;
corrected_distance = distance - drift_correction;


% Time array for X-axis in plotting
sample_interval = 40e-6; % Sampling interval (40 microseconds)
total_samples = num_cols * num_frm_con;
time = (1:total_samples) * sample_interval;

% Plot Relative Distance
figure('Name', 'Relative Distance');
plot(time, corrected_distance);
title('Relative Distance');
xlabel('time [s]');
ylabel('Distance [mm]');

% Calculate Cumulative (Absolute) Distance
% Using moving median filter to smooth the distance signal
smoothed_distance = movmedian(corrected_distance, 5000);


% Calculating cumulative distance over time
cumulative_distance = cumsum(smoothed_distance);

% Plot Absolute Distance
figure('Name', 'Cumulative distance');
plot(time, cumulative_distance);
title('Cumulative Distance');
xlabel('time [s]');
ylabel('Distance [mm]');

checksum = sum(cumulative_distance(:));
disp(['Checksum after Calculating cumulative distance: ', num2str(checksum)]);


% Save the cumulative distance data to a .mat file
save(fullfile(data_files_folder,[case_file,'_cumulative_distance.mat']), 'cumulative_distance');

%% below has not been refactored yet 1/4/2024. 
% Question what is the purpose of the following code?
% Looks like publication quality image

figure(11)
% Create axes
axes1 = axes('Parent',figure(11),...
    'Position',[0.623353293413173 0.167330677290833 0.449101796407186 0.697211155378486]);
hold(axes1,'on');
alpha 0.2

plot(time,totaldisr(:,:),'--o','Color',[0.77,0.72,0.97],'LineWidth',4,'MarkerSize',10,'MarkerFaceColor',[0.47,0.36,0.94],'MarkerEdgeColor',[0.47,0.36,0.94]);
%figure(11), legend({'Observed','Predicted'},'LineWidth',2);
figure(11), title('Doppler-tracked Absolute Distance','FontSize',40);
figure(11), xlabel('time [s]','FontSize',36);
figure(11), ylabel('Distance [mm]','FontSize',36);
xlim([0,25])
% ylim([-5,55])
set(gca,'FontSize',36);
set(gcf, 'Units', 'normalized', 'outerposition', [0, 0, 1, 0.7],'Color',[1 1 1],'OuterPosition',[0 0 1 1])
set(get(gca,'YLabel'),'Rotation',0)
ylh = get(gca,'ylabel');
gyl = get(ylh);                
ylp = get(ylh, 'Position');
ylp(1) = -5;
set(ylh, 'Rotation',90, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','center')
ax = gca;
fig_loc = ax.Position;
ax.Position(1) = 0.1802;
pbaspect([1 1 1]);
% set(axes1,'LineWidth',6,'TickLength',[0.025 0.04],'XTick',...
%     [0 5 10 15 20 25],'FontSize',36);
%% below has not been refactored yet 1/4/2024. 
% Question what is the purpose of the following code?
% likely to be an older version of distance calculation snippet above

%Distance calibration line
w_length=0.0013; %mm %1.3um
n=1.3;
disc=(diff_phase_con_med(num_rows_cal,:,:).*w_length)./(4*pi*n);
disc=disc-mean(disc(:,1:10000));
time=(1:num_cols*num_frm_con)*40*(10^-6);
figure('Name','Relative distance c'),plot(time,disc(:,:));
title('Relative Distance')
xlabel('time [s]') 
ylabel('Distance [mm]')
sum=0;
totaldisc=zeros(1,num_cols*num_frm_con);
for u=1:num_cols*num_frm_con
    sum=sum+disc(1,u);
    totaldisc(:,u)=sum;
end
figure('Name','Absolute Distance c'),plot(time,totaldisc(:,:));
title('Absolute Distance')
xlabel('time [s]') 
ylabel('Distance [mm]')


% %Tiff file
% filename = 'phasedifference.tiff';
% for k=1:num_frm_con
%     B=diffphconM1(:, :, k);
%     %figure,mesh(B)
%     B=uint8(B.*(255/(2*pi)));
%     B=ind2rgb(B,redblue);
%     if k==1
%         imwrite(B,filename);
%     else
%         imwrite(B,filename,'WriteMode','append');
%     end
%     
% end


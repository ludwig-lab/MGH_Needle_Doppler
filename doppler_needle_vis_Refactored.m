%This code was developed to use pre-stored phase differences measured from 
%a needle probe and create surface-referenced maps of intensity, phase 
%retardation and optic axis as the needle moves.  

%V1.0, Danielle J. Harper and Benjamin J. Vakoc, 2023
%Further information can be found in this paper:
   % https://doi.org/10.48550/arXiv.2305.14390
clc;
close all;
clear all;

case_file = 'example_salmon'; 

% Switch between cases
switch case_file
    case 'example_salmon' % Example data set       
        data_files_folder = 'E:\My Drive\4_UW-Madison\2_4_MGH+FDA+UW OCT\needle_code_share\Example_data';
        ii_file = '[p.needle_2][s.salmon][06-16-2021_14-15-13]ii.mgh';
        oa_avg_file = 'OAavg_salmon.tif';
        rr_avg_file = 'rravg_salmon.tif';
        ii_avg_file = 'iiavg_salmon.tif';
        machine_ID = 'MGH';

    case 'UW_salmon' % Salmon data
        data_files_folder = 'H:\OFDIData\user.Ricardo\[p.231220_PS_Needle_Probe]\[p.231220_PS_Needle_Probe][s.Salmon_Probe_04_Test3_Dist1.2cm][12-20-2023_15-52-35]';
        ii_file = '[p.231220_PS_Needle_Probe][s.Salmon_Probe_04_Test3_Dist1.2cm][12-20-2023_15-52-35].structure8.mgh';
        oa_avg_file = 'AVG_Reslice_oa1.tif'; 
        rr_avg_file = 'AVG_Reslice_rr.tif';
        ii_avg_file = 'AVG_Reslice_ii.tif';
        machine_ID = 'SPARC';

    otherwise
        error('Invalid case name. Please specify `example_salmon` or `UW_salmon`.');
end

% File loading operations
ii_file_id = fopen(fullfile(data_files_folder, ii_file), 'rb'); % Open intensity file
oa_avg = imread(fullfile(data_files_folder, oa_avg_file));
rr_avg = imread(fullfile(data_files_folder, rr_avg_file));
ii_avg = imread(fullfile(data_files_folder, ii_avg_file));
cumulative_distance_file = fullfile(data_files_folder, [case_file, '_cumulative_distance.mat']);
load(cumulative_distance_file); % Load the cumulative distance data

num_rows = size(ii_avg,1);
num_cols = size(ii_avg,2);

%%
z = cumulative_distance(1:511:end)*600;
pixel = round(z/max(z)*1399)+1;
pixel(pixel<1)=1;
%needle visualization
needle_visualization_matrix = [-10 0; 10 0; 10 1425; -10 1385] ; 

fseek(ii_file_id,1024*1024,'bof'); % skipping the meta data

angle_range = 88:0.05:92;
v = [];
clear rr_Int;
oa_Int = uint8(255*(rand(3600,200,3)*0.5+0.5));
rr_Int = uint8(255*(rand(3600,200)*0.5+0.5));
ii_Int = uint8(255*(zeros(3600,200)));

height = zeros(3600,1,1);
writemode = 'overwrite';

 
writerObj = VideoWriter('salmon_rr.avi','Uncompressed AVI');
 writerObj.FrameRate=50;
 open(writerObj);
%%
sampling_interval = 40e-6; % 40 microseconds, used earlier in your script

% New constants for processing
scaling_factor = 600; % Scaling factor used to scale cumulative_distance for normalization
max_pixel_value = 1399; % Maximum pixel value, used for normalizing pixel intensities
needle_visualization_matrix = [-10 0; 10 0; 10 1425; -10 1385]; % Coordinates for needle visualization
angle_range = 88:0.05:92; % Range of angles for image processing or analysis
image_dimensions = [3600, 200]; % Dimensions (rows, columns) for the image arrays
frame_rate = 50; % Frame rate for the output video
video_name = 'salmon_rr.avi'; % Name of the output video file

% Calculations for image processing
num_cols_phase_diff_arr = floor(num_cols/2)-1; % Calculate number of columns for phase difference array
z = cumulative_distance(1:num_cols_phase_diff_arr:end) * scaling_factor; % Scale the cumulative_distance by the scaling factor
pixel = round(z / max(z) * max_pixel_value) + 1; % Normalize and scale `z`, then round and offset by 1 to avoid zero-indexing?
pixel(pixel < 1) = 1; % Ensure pixel values do not fall below 1, #Question: what does negative mean here?


% File operations
fseek(ii_file_id, num_rows * num_cols, 'bof'); % Move file position indicator to the start of the desired data

% Initialization of variables and arrays
v = []; % Initialize an empty array for storing processed values
clear rr_Int; % Clear any existing data in rr_Int
oa_Int = uint8(255 * (rand(image_dimensions(1), image_dimensions(2), 3) * 0.5 + 0.5)); % Initialize oa_Int with random values for optical axis visualization
rr_Int = uint8(255 * (rand(image_dimensions(1), image_dimensions(2)) * 0.5 + 0.5)); % Initialize rr_Int with random values for retardation visualization
ii_Int = uint8(zeros(image_dimensions(1), image_dimensions(2))); % Initialize ii_Int with zeros for intensity visualization

height = zeros(image_dimensions(1), 1, 1); % Initialize height array for storing height information
writemode = 'overwrite'; % Set file write mode to overwrite for video output

% Setup for writing the video
writerObj = VideoWriter(video_name, 'Uncompressed AVI'); % Create a VideoWriter object for the output video
writerObj.FrameRate = frame_rate; % Set the frame rate for the video
open(writerObj); % Open the video file for writing

 
 %% 
% Need to change from the default renderer to zbuffer to get it to work right.
% openGL doesn't work and Painters is way too slow.
set(gcf, 'renderer', 'zbuffer');


for i=1:700
    
    fseek(ii_file_id,1024*1024 + 1024*1024*(i-1+200),'bof');
    img = fread(ii_file_id,[1024 1024],'uint8=>single').';
    %imgrr = fread(fidrr,[1024 1024],'uint8=>single').';
    
    for j=1:8
        
        simg = img(500:600,(1:128) + (j-1)*128);
        R = radon(edge(simg),angle_range);
        [~,mi] = max(std(R));
        v = [v angle_range(mi)];
        
    end
   % temp = cumsum(tand((v-90))*128);
    temp = pixel(i+200);
    pos(i) = round(temp(end));
 %    pos(i) = z(i)+1;
    %rr(i) = mean(mean(imgrr(500:600,:)));
%     figure(1);
%     plot(temp);
%     drawnow;
    pos = (pos + 1401);
    
    figure(2);
    height(pos(i):(pos(i)+554)) = height(pos(i):(pos(i)+554)) + 1;
    oa_Int(pos(i):(pos(i)+554),2:end,:) = oa_Int(pos(i):(pos(i)+554),1:(end-1),:);
    oa_Int(pos(i):(pos(i)+554),1,:) = oa_avg(470:1024,i+200,:);
    
    rr_Int(pos(i):(pos(i)+554),2:end,:) = rr_Int(pos(i):(pos(i)+554),1:(end-1));
    rr_Int(pos(i):(pos(i)+554),1,:) = rr_avg(470:1024,i+200);
    
    ii_Int(pos(i):(pos(i)+554),2:end,:) = ii_Int(pos(i):(pos(i)+554),1:(end-1));
    ii_Int(pos(i):(pos(i)+554),1,:) = ii_avg(470:1024,i+200);
    
    
    oa_sInt = oa_Int((pos(i)-1400):(pos(i)+400),:,:);
    oa_sInt(1400:1410,:,:) = 0;
    
    rr_sInt = rr_Int((pos(i)-1400):(pos(i)+400),:);
    rr_sInt(1400:1410,:) = 0;
    
    ii_sInt = ii_Int((pos(i)-1400):(pos(i)+400),:);
    ii_sInt(1400:1410,:) = 0;
    
    height_s = height((pos(i)-1400):(pos(i)+400));
    height_s(height_s>199)=199;
    %if (mod(i,10)==0)
        oa_sIntB = zeros(length(height_s),200,3,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                oa_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
                oa_sIntB(j,1:200,2) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),2)),linspace(1,sqrt(height_s(j)),200).^2));
                oa_sIntB(j,1:200,3) = uint8(interp1(1:height_s(j),single(oa_sInt(j,1:height_s(j),3)),linspace(1,sqrt(height_s(j)),200).^2));
            else
                oa_sIntB(j,1:200,1:3) = uint8(0);
            end
        end
        
        rr_sIntB = zeros(length(height_s),200,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                rr_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(rr_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
               
            else
                rr_sIntB(j,1:200) = uint8(0);
                
            end
        end
        
                ii_sIntB = zeros(length(height_s),200,'uint8');
        for j=1:length(height_s)
            if (height_s(j)>1)
                ii_sIntB(j,1:200,1) = uint8(interp1(1:height_s(j),single(ii_sInt(j,1:height_s(j),1)),linspace(1,sqrt(height_s(j)),200).^2));
               
            else
                ii_sIntB(j,1:200) = uint8(0);
                
            end
        end
        
     cmap=colormap(bone(256));
     cmap = cmap(:,[2 1 3]);
    %choose which plot to make here - ii, rr or oa
       % imshow(imrotate(imfuse(flip(ii_sIntB,2),ii_sIntB,'Montage'),0));
        imshow(imrotate(imfuse(flip(rr_sIntB,2),rr_sIntB,'Montage'),0),cmap);
       % imshow(imrotate(imfuse(flip(oa_sIntB,2),oa_sIntB,'Montage'),0));
        %imwrite(imrotate(imfuse(flip(oa_sIntB,2),oa_sIntB,'Montage'),90),'out.tif','WriteMode',writemode);writemode = 'append';
        
        patch([0 0 400 400],[1400 1800 1800 1400],'y','FaceAlpha',0.1)
        patch(200+needle_visualization_matrix(:,1),needle_visualization_matrix(:,2),[0.5,0.5,0.5])
       %   patch(200+needle_visualization_matrix(:,1),pixel(i+200)-117+needle_visualization_matrix(:,2),'c')
        
        drawnow;
        i
        
        thisFrame = getframe(gca);
	% Write this frame out to a new video file.
  	writeVideo(writerObj, thisFrame);
	%myMovie(frameIndex) = thisFrame;
    %end
    
    
    
end
close(writerObj)
fclose(ii_file_id);
        

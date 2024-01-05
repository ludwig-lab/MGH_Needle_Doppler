function diff_ph_con = concatenate_frames(diff_ph1_sliced)
    % concatenate_frames Concatenates pairs of frames in a 3D array.
    %
    % This function concatenates every two frames in a 3D array, handling cases 
    % where the number of frames might be odd. In the case of an odd number of 
    % frames, the last frame is added as is, without concatenation.
    %
    % Inputs:
    %   diff_ph1_sliced: A 3D array of sliced phase differences.
    %
    % Outputs:
    %   diff_ph_con: A 3D array containing the concatenated frames.

    % Initialize parameters for concatenation
    num_rows = size(diff_ph1_sliced,1); %   num_rows: The number of rows in the diff_ph1_sliced array.
    num_frm = size(diff_ph1_sliced,3); %   num_frm: The total number of frames in the diff_ph1_sliced array.
    frameStep = 2; % Step size for iterating through frames
    halfFrameWidth = size(diff_ph1_sliced, 2); % Half the width of the frame
    concatenationWidth = halfFrameWidth * 2; % Total width for concatenation

    % Initialize the concatenated phase difference array
    diff_ph_con_ = zeros(num_rows, concatenationWidth, num_frm);

    % Iterate through the frames with the specified step size
    for l = 1:frameStep:num_frm
        if l + 1 > num_frm
            % Handle the case where the frame number exceeds the limit
            diff_ph_con_(:, 1:halfFrameWidth, l) = diff_ph1_sliced(:, :, l);
            break;
        end
        % Concatenate two frames together
        diff_ph_con_(:, :, l) = horzcat(diff_ph1_sliced(:, :, l), diff_ph1_sliced(:, :, l + 1));
    end

    % Slice the concatenated array to remove frames with only zeros
    diff_ph_con = diff_ph_con_(:, :, 1:frameStep:end);
end

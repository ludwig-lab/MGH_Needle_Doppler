function diff_ph1 = calculate_phase_difference(ph1, col_distance)
    % calculate_phase_difference Calculates phase differences in OCT data.
    %
    % This function computes the phase differences in Optical Coherence
    % Tomography (OCT) data based on a specified column distance. It is designed
    % to handle different types of OCT data represented by the column distance.
    % For instance, a column distance of 2 is used for 'MGH' type data where
    % the phase difference is calculated by skipping one column, while a
    % column distance of 1 is used for 'SPARC' type data for calculating the
    % phase difference between adjacent columns.
    %
    % Inputs:
    %   ph1: A 3D matrix representing the OCT phase data. The dimensions of
    %        this matrix are expected to be [rows, columns, frames].
    %   col_distance: An integer value representing the column distance for
    %                 calculating the phase difference. A value of 2 indicates
    %                 skipping one column (for 'MGH'), and 1 indicates
    %                 adjacent columns (for 'SPARC').
    %
    % Outputs:
    %   diff_ph1: A 3D matrix of the same size as 'ph1'. This matrix contains
    %             the calculated phase differences. The dimensions of the matrix
    %             are adjusted based on the column distance parameter.
    %
    % Example:
    %   diff_ph1_mgh = calculate_phase_difference(ph1_data, 2);
    %   diff_ph1_sparc = calculate_phase_difference(ph1_data, 1);

    % Calculate the size of the output array based on col_distance
    [rows, cols, f] = size(ph1);

    diff_ph1 = zeros(rows, floor(cols/col_distance) - 1, f);


    % Computing phase differences
    for j = 1:f
        % Iterate over frames
        for i = 1:col_distance:cols - col_distance
            % Calculate the phase difference based on column distance
            diff = ph1(:, i, j) - ph1(:, i + col_distance, j);

            % Store the computed phase difference
            if col_distance == 2
                % Store in corresponding index for MGH
                diff_ph1(:, (i + 1)/2, j) = diff;
            else
                % Store in corresponding index for SPARC
                diff_ph1(:, i, j) = diff;
            end
        end
    end
end

function crc = crc32(data)
    % crc32 Computes the CRC-32 checksum for various data types.
    % 
    % Syntax:
    %   crc = crc32(data)
    %
    % Description:
    %   This function computes the CRC-32 checksum of the given data. It accepts
    %   data in various formats including both floating-point and integer arrays.
    %   The function converts non-uint8 data types into a byte stream using
    %   typecasting, which ensures that the binary representation of the data
    %   is accurately processed for the checksum calculation.
    %
    % Inputs:
    %   data - The input data for which the CRC-32 checksum is to be computed.
    %          This can be an array of any data type (e.g., uint8, single, double).
    %
    % Outputs:
    %   crc  - The computed CRC-32 checksum, returned as a uint32 value.
    %
    % Examples:
    %   % Example 1: Compute CRC-32 for a vector of uint8 data
    %   crc1 = crc32(uint8([123, 224, 56, 78]));
    %
    %   % Example 2: Compute CRC-32 for a vector of floating-point numbers
    %   crc2 = crc32([1.23, 2.34, 3.45]);
    %
    %   % Example 3: Compute CRC-32 for a matrix of mixed data types
    %   crc3 = crc32(single([1.2, 3.4; 5.6, 7.8]));
    %
    % Note:
    %   For non-uint8 data types, the function uses MATLAB's 'typecast'
    %   function to convert the data into a byte stream. This conversion
    %   process ensures that all bits of the data's binary representation
    %   are considered in the CRC-32 calculation.
    %
    %   The CRC-32 algorithm uses a standard polynomial represented in
    %   hexadecimal as 'EDB88320'. The algorithm processes each byte of
    %   the data, updating the checksum accordingly.

    % Initialize variables
    crc  = uint32(hex2dec('FFFFFFFF'));
    poly = uint32(hex2dec('EDB88320'));

    % Check data type and convert to byte stream if necessary
    if ~isa(data, 'uint8')
        data = typecast(data(:), 'uint8');
    end

    % Compute CRC-32 value
    for i = 1:length(data)
        crc = bitxor(crc, uint32(data(i)));
        for j = 1:8
            mask = bitcmp(bitand(crc, uint32(1)));
            if mask == intmax('uint32'), mask = 0; else mask = mask + 1; end
            crc = bitxor(bitshift(crc, -1), bitand(poly, mask));
        end
    end
    crc = bitcmp(crc);
end
function [ qamSymbols ] = mapper( tx_bits,b,N_frame)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%   tx_bits have n rows of 4 bits
%   1. we first find the number of rows in the input matrix
%   2. we loop each row to read the bits value in row
    reshapeMatrix = reshape(tx_bits,[b,N_frame*2]);
    reshapeMatrix = reshapeMatrix';
    row = size(reshapeMatrix);
    %init symbols array to keep all symbol after convert bits frame to symbol.
    qamSymbols = [];
    for rowCounter = 1:row
        frame = reshapeMatrix(rowCounter,:);

        % convert x to a string array
        bitToString = num2str(frame);

        % remove the sp aces in the string, so it becomes '0000'
        bitToString(isspace(bitToString)) = '';
        %build mapper
        switch bitToString
            case '0000'
                encodedSymbol =  -3-3*j;
            case '0001'
                encodedSymbol =  -3-j;
            case '0010'
                encodedSymbol =  -3+3*j;
            case '0011'
                encodedSymbol =  -3+j;
            case '0100'
                encodedSymbol =  -1-3*j;
            case '0101'
                encodedSymbol =  -1-j;
            case '0110'
                encodedSymbol =  -1+3*j;
            case '0111'
                encodedSymbol =  -1+j;
            case '1000'
                encodedSymbol =  3-3*j;
            case '1001'
                encodedSymbol =  3-j;
            case '1010'
                encodedSymbol =  3+3*j;
            case '1011'
                encodedSymbol =  3+j;
            case '1100'
                encodedSymbol =  1-3*j;
            case '1101'
                encodedSymbol =  1-j;
            case '1110'
                encodedSymbol =  1+3*j;
            case '1111'
                encodedSymbol =  1+j;
            otherwise
        end
        % merge symbol to the symbols array
        qamSymbols(rowCounter)= encodedSymbol;
    end
    qamSymbols = reshape(qamSymbols,[rowCounter,1]);
end
function [ softDecodedBits ] = softdemapper( rx_sym)
%DEMAPPER Summary of this function goes here
%   Detailed explanation goes here
    % real number to binary
    softDecodedBits = [];
    %rx_sym
    % we loop each row to read the encoded symbol in row
    row = size(rx_sym,1);
    for symbolRow = 1:row
        symbol = rx_sym(symbolRow,:);
        softBit = softbit(symbol);
        softDecodedBits = [softDecodedBits softBit];
    end
    softDecodedBits = softDecodedBits.';
end

function [soft_qam_symbol] = softbit(realNumber)
    Yre = real(realNumber); % real part
    Yim = imag(realNumber); % imaginary part
    b0 = 0;
    b1 = 0;
    b2 = 0;
    b3 = 0;
    % b0 = 2(Yre + 1) if Yre<-2
    % b0 = Yre if -2 <= Yre < 2
    % b0 = 2(Yre - 1) Yre > 2
    
    % b1 = -|Yre|+2 for all Yre
    
    % b2 = 2(Yim + 1) if Yim<-2
    % b2 = Yim if -2 <= Yim < 2
    % b2 = 2(Yim - 1) Yim > 2
 
    % b3 = -|Yim|+2 for all Yim
    
    if Yre < -2
        b0 = 2*(Yre + 1);
    elseif Yre >= -2 & Yre < 2
        b0 = Yre;
    elseif Yre >=2
        b0 = 2*(Yre - 1);
    end
    
    b1 = -abs(Yre)+2;
    
    if Yim < -2
        b2 = 2*(Yim + 1);
    elseif Yim >= -2 & Yim < 2
        b2 = Yim;
    elseif Yim >=2
        b2 = 2*(Yim - 1);
    end
    
    b3 = -abs(Yim)+2;
    
    bit0 = (sign(b0)+1)/2;
    bit1 = (sign(b1)+1)/2;
    bit2 = (sign(b2)+1)/2;
    bit3 = (sign(b3)+1)/2;
    soft_qam_symbol = [bit0 bit1 bit2 bit3];
end


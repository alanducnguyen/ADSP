function [ hardDecodedBits ] = harddemapper( rx_sym)
%DEMAPPER Summary of this function goes here
%   Detailed explanation goes here
    % real number to binary
    hardDecodedBits = [];
    %rx_sym
    % we loop each row to read the encoded symbol in row
    row = size(rx_sym,1);
    for symbolRow = 1:row
        symbol = rx_sym(symbolRow,:);
        hardBit = hardbit(symbol);
        hardDecodedBits = cat(1,hardDecodedBits,hardBit);
    end
end

function [qam_symbol] = hardbit(realNumber)
    Yre = real(realNumber); % real part
    Yim = imag(realNumber); % imaginary part
    % range Yre <=-2, then Yim {[<=-2],[-2,0],[0,2],[>2]}
    if Yre <= -2
        if Yim <=-2
            qam_symbol = [0 0 0 0];
        end
        
        if Yim > -2 & Yim <=0
            qam_symbol = [0 0 0 1];
        end
        
        if Yim >0 & Yim <=2
            qam_symbol = [0 0 1 1];
        end
        
        if Yim >2
            qam_symbol = [0 0 1 0];
        end
    end
    % range Yre >-2 and Yre <=0, then Yim {[<=-2],[-2,0],[0,2],[>2]}
    if Yre > -2 & Yre <=0 
        if Yim <=-2
            qam_symbol = [0 1 0 0];
        end
        
        if Yim > -2 & Yim <=0
            qam_symbol = [0 1 0 1];
        end
        
        if Yim >0 & Yim <=2
            qam_symbol = [0 1 1 1];
        end
        
        if Yim >2
            qam_symbol = [0 1 1 0];
        end
    end
    % range Yre >0 and Yre <=2, then Yim {[<=-2],[-2,0],[0,2],[>2]}
    if Yre > 0 & Yre <=2 
        if Yim <=-2
            qam_symbol = [1 1 0 0];
        end
        
        if Yim > -2 & Yim <=0
            qam_symbol = [1 1 0 1];
        end
        
        if Yim >0 & Yim <=2
            qam_symbol = [1 1 1 1];
        end
        
        if Yim >2
            qam_symbol = [0 1 1 0];
        end
    end
    % range Yre >2 then Yim {[<=-2],[-2,0],[0,2],[>2]}
    if Yre > 2
        if Yim <=-2
            qam_symbol = [1 0 0 0];
        end
        
        if Yim > -2 & Yim <=0
            qam_symbol = [1 0 0 1];
        end
        
        if Yim >0 & Yim <=2
            qam_symbol = [1 0 1 1];
        end
        
        if Yim >2
            qam_symbol = [1 0 1 0];
        end
    end
end
function [ outputBits ] = myconvcode( inputData )
    %MYCONVCODE Summary of this function goes here
    %   Detailed explanation goes here
    % s0 s1 s2
    registerMatrix = zeros(1,3);
    outputBits = [];
    for i=1:length(inputData)
        b = inputData(i);

        % shift stageMatrix to the left 1 stage, 
        % then input bit b = s0, s1 = prev s0, s2 = prev s1
        registerMatrix = [b registerMatrix(1:end-1)];

        %get current state from register matrix
        currentState = registerMatrix(2:end);

        outputResult = ComputeStateTransition(b,currentState);
        outputBits = [outputBits outputResult];
    end
end



% function output = ComputeStateTransition(b,currentState)
%     s0 = b;
%     s1 = currentState(1,1);
%     s2 = currentState(1,2);
%     
%     bit0 = mod(s0+s2,2);
%     bit1 = mod(s0 + s1 + s2,2);
%     output = [bit0 bit1];
% end

function output = ComputeStateTransition(b,currentState)
if b == 0
if isequal([0 0], currentState)
    output =  [0 0];
end
if isequal([0 1], currentState)
    output =  [1 1];
end
if isequal([1 0], currentState)
    output =  [1 0];
end
if isequal([1 1], currentState)
    output =  [0 1];
end
end
if b == 1
if isequal([0 0], currentState)
    output =  [1 1];
end
if isequal([0 1], currentState)
    output =  [0 0];
end
if isequal([1 0], currentState)
    output =  [0 1];
end
if isequal([1 1], currentState)
    output =  [1 0];
end
end

end


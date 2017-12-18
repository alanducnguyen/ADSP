function [ convData ] = convencode(inputData)
    %CONVCODE Summary of this function goes here
    %   n = k + r
    %   k = information bits
    %   r = extra bits
    %   R = k/n, R high then the result will be more exactly.
    %   N = constraint length (memory of N+1 stages)
    %   k = 1 mean that, each time we send 1 bit.
    %   N = 2 mean that half-rate constraint length.
    %   g0(x) = 1 + x^2
    %   g1(x) = 1 + x + x^2
    
    %  s0 s1 s2
    convData = [];
    registerMatrix = zeros(1,3);
    for i=1:length(inputData)
        b = inputData(i);
        % shift stageMatrix to the left 1 stage, 
        % then input bit b = s0, s1 = prev s0, s2 = prev s1
        registerMatrix = [b registerMatrix(1:end-1)];
        %get current state from register matrix
        currentState = registerMatrix(2:end);
        outputResult = ComputeStateTransition(b,currentState);
        convData = [convData outputResult];
    end
end

function output = ComputeStateTransition(b,currentState)
    s0 = b;
    s1 = currentState(1,1);
    s2 = currentState(1,2);
    % modulo 2 adder
    % bit0 = modulo 2 adder of s0 and s2
    bit0 = mod(s0+s2,2);
    
    % bit1= modulo 2 adder of s0 ,s1 and s2
    modStep1 = mod(s0+s1,2);
    bit1 = mod(modStep1+s2,2);
    output = [bit0 bit1];
end


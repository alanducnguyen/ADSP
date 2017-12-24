function [ decoded_bits ] = myviterbi( input_bits )
% input_bits = [
%                 0 0
%                 1 1
%                 1 1
%                 0 0
%                 0 1
%                 1 0
%                 0 1
%                 1 1
%                 1 1
%                 1 0
%                 0 0
%                 0 0
%                 1 1
%                 0 0
%                 1 1
%                 1 0
%                 1 1
%              ];
%MYVITERBI Summary of this function goes here
%   Detailed explanation goes here
% INPUT:    00 11 10 00 | 01 10 01 11 | 11 10 00 10 | 11 00 11 10 | 11         
% RECEIVED: 00 11 11 00 | 01 10 01 11 | 11 10 00 00 | 11 00 11 10 | 11 
% WRONG:          x                          x
% STATE: table next state
% (0)00 | 0 0 2 
% (1)01 | 1 0 2
% (2)10 | 2 1 3
% (3)11 | 3 1 3
% ns_table=[0 0 2; 1 0 2; 2 1 3; 3 1 3]; %NEXT STATE array
% END STATE
% --------------------
% POSSIBLE OUTPUT: table output
% output_table=[00 00 11; 01 11 00; 10 10 01; 11 01 10]; %OUTPUT array
% --------------------
% steps:
% foreach rxBits length(max steps)
% as the register always init with 0 0 0, 
% then the first bit we can receive
% alway 00 or 11: ref table output, 
% current state = 00 then next state posible is: 00 or 10, ref: table next
% steps
% 
    N = length(input_bits);
    % keep the weight based on each steps
    distance(1:4, 1:N+1)=Inf;
    % keep the parent state node history
    state_history(1:4, 1:N+1)=Inf;

    % loop N times to compute to find the winning path.
    for st = 1:1:N+1
        % if st = 1
        if st == 1
            %Initialization of distance for start point = 0
            distance(1,1)=1;
            % set the previous position at 1 (state  = 00)
            previous_pos=[1];
            % set state history = 0, this variable will keep the parent of each
            % posible next state
            state_history(1,1) = 0;
        else
            % loop each previous(parent) positions, then get the next transision pos(next state)
            temp_state=[];
            % with each posible root node, then we can retrieve the posible
            % state, assume input = 0 or = 1
            for parent_pos=1:1:length(previous_pos)
                temp_pos = previous_pos(1,parent_pos);
                % get next state if the input bit = 0
                next_state_0 = getNextState(0,temp_pos);
                % get next state if the input bit = 1
                next_state_1 = getNextState(1,temp_pos);
                % get output if the input bit = 0
                outputO=getOutputBit(0, temp_pos);
                % get output if the input bit = 1
                output1=getOutputBit(1, temp_pos);
                % get the pair bits symbol at step st-1
                step_bit = input_bits(st-1,:);
                %calculate bit error if bit = 0
                weight_at_0 = sum(abs(step_bit - outputO));
                %calculate bit error if bit = 1
                weight_at_1 = sum(abs(step_bit - output1));
                %retrieve the latest weight in the recorded distances matrix
                parent_weight = distance(temp_pos,st-1);

                if parent_weight == 'Inf'
                    % there is nothing to do, if parent weight = Inf
                else
                    % in case weight of position has not been recorded before
                    % then we record it as first value
                    if distance(next_state_0,st) == 'Inf'
                        distance(next_state_0,st) = weight_at_0 + parent_weight;
                    else
                        % case 0: if the value has been set, then we have to determine
                        % if the value is minimum yet? if there is any smaller
                        % value, then replace this value as weight of tracking
                        % position
                        current_0_value = distance(next_state_0,st);
                        new_0_value = weight_at_0 + parent_weight;
                        if new_0_value<current_0_value
                            distance(next_state_0,st) = new_0_value;
                            state_history(next_state_0,st) = temp_pos;
                        end
                    end

                    % in case weight of position has not been recorded before
                    % then we record it as first value
                    if distance(next_state_0,st) == 'Inf'
                        distance(next_state_0,st) = weight_at_0 + parent_weight;
                    else
                        % case 1: if the value has been set, then we have to determine
                        % if the value is minimum yet? if there is any smaller
                        % value, then replace this value as weight of tracking
                        % position
                        current_1_value = distance(next_state_1,st);
                        new_1_value = weight_at_1 + parent_weight;
                        if new_1_value<current_1_value
                            distance(next_state_1,st) = weight_at_1 + parent_weight;
                            state_history(next_state_1,st) = temp_pos;
                        end
                    end
                end
                % keep the temp_state = every posible state, then it will be
                % the parents node for next steps
                temp_state = [temp_state next_state_0 next_state_1];
            end
            % END FOR previous_pos
            % reset previous_pos = all new posible state, then this value will
            % become parent node for next steps.
            previous_pos = unique(temp_state);

        end

        distance_flip = fliplr(distance);
        % END checking st
    end
%     distance_flip =[
%      4     3     3     2     3     2     2     1     0
%      4     4     3     3     2     3     1   Inf   Inf
%      4     3     3     3     3     2     2     1   Inf
%      3     4     3     3     2     2     3   Inf   Inf
%         ];
%     distance_flip
%     smallest_weight = min(distance_flip(:, 8))
%     
%     smallest_position = find( distance_flip(:, st)==smallest_weight )
%     
%     return
    
    % winning path sequences
    winning_path = ones(1, N); %STATE SEQUENCE array
    for st = 1:1:N+1
        %find the smallest value in recorded distance
        smallest_weight = min(distance_flip(:, st));
        smallest_position = find( distance_flip(:, st)==smallest_weight );
        if(length(smallest_position)>1 && st>=2)
            previous = winning_path(st-1);
            for option=1:1:length(smallest_position)
                pos_option = smallest_position(option);
                if (previous == 1) && (pos_option==1 || pos_option==2)
                    winning_path(st)=pos_option;
                end
                if (previous == 2) && (pos_option==3 || pos_option==4)
                    winning_path(st)=pos_option;
                end
                if (previous == 3) && (pos_option==1 || pos_option==2)
                    winning_path(st)=pos_option;
                end
                if (previous == 4) && (pos_option==3 || pos_option==4)
                    winning_path(st)=pos_option;
                end
            end
        else
            winning_path(st)=smallest_position(1);
        end
    end
    
    % revert the winning path to get the correct order.
    winning_path = fliplr(winning_path);
    decoded_bits = [];
    for state = 1:1:length(winning_path)-1
        start_state = winning_path(1,state);
        next_state = winning_path(1,state+1);
        bit = get_bit_trace_back(start_state,next_state);
        decoded_bits = [decoded_bits bit];
    end
    % END loop N times to compute to find the winning path.
end

function [nextState] = getNextState(inputBit,currentState)
% column 1 = current state
% column 2 = output if input bit = 0
% column 3 = output if input bit = 1
% (1)00 | 1 1 3 
% (2)01 | 2 1 3
% (3)10 | 3 2 4
% (4)11 | 4 2 4
    nextState = Inf;
    if inputBit == 0
        if isequal(currentState,1)
            nextState = 1;
        end

        if isequal(currentState,2)
            nextState = 1;
        end

        if isequal(currentState,3)
            nextState = 2;
        end

        if isequal(currentState,4)
            nextState = 2;
        end
    end

    if inputBit == 1
        if isequal(currentState,1)
            nextState = 3;
        end

        if isequal(currentState,2)
            nextState = 3;
        end

        if isequal(currentState,3)
            nextState = 4;
        end

        if isequal(currentState,4)
            nextState = 4;
        end
    end
end

function [outputBit] = getOutputBit(inputBit, currentState)
% column 1 = current state
% column 2 = output if input bit = 0
% column 3 = output if input bit = 1
% op_table=[00 00 11; 01 11 00; 10 10 01; 11 01 10]; %OUTPUT array
    outputBit = Inf;
    if inputBit == 0
        if isequal(currentState,1)
            outputBit = [0 0];
        end

        if isequal(currentState,2)
            outputBit = [1 1];
        end

        if isequal(currentState,3)
            outputBit = [1 0];
        end

        if isequal(currentState,4)
            outputBit = [0 1];
        end
    end

    if inputBit == 1
        if isequal(currentState,1)
            outputBit = [1 1];
        end

        if isequal(currentState,2)
            outputBit = [0 0];
        end

        if isequal(currentState,3)
            outputBit = [0 1];
        end

        if isequal(currentState,4)
            outputBit = [1 0];
        end
    end
end

function [bit_trace_back] = get_bit_trace_back(state,next_state)
    % here are all possible transition state, vertical then up
    %        |   00(1)   01(2)   10(3)   11(4)
    % -------------------------------------------
    % 00(1)  |   0       Inf     1       Inf     
    % 01(2)  |   0       Inf     1       Inf
    % 10(3)  |   Inf     0       Inf     1   
    % 11(4)  |   Inf     0       Inf     1
    %--------------------------------------------
    %
    % transition_status= [0 Inf 1 Inf; 0 Inf 1 Inf; Inf 0 Inf 1; Inf 0 Inf 1];
    transition_status= [0 Inf 1 Inf; 0 Inf 1 Inf; Inf 0 Inf 1; Inf 0 Inf 1];
    %state
    %next_state
    bit_trace_back = transition_status(state,next_state);
    if bit_trace_back == Inf
        bit_trace_back=round(rand);
    end
end

function [State State_change1 State_change0] = dictionary()
State = zeros(16,2);
for i = 1:16
    if mod(i,4) == 1
        State(i,1) = 3;
    elseif mod(i,4) == 2
           State(i,1) = 1;
    elseif mod(i,4) == 3
           State(i,1) = -1;
    else 
           State(i,1) = -3;    
    end
    if floor((i-1)/4) == 0
        State(i,2) = 3;
    elseif floor((i-1)/4) == 1
        State(i,2) = 1;
    elseif floor((i-1)/4) == 2
        State(i,2) = -1;
    else
        State(i,2) = -3;
    end   
end
for i = 1:16
    for j = 1:2
        if State(i,j) == 3
            State_change{i,j} = [1,1];
        elseif State(i,j) == 1
            State_change{i,j} = [1,-1];
        elseif State(i,j) == -1
            State_change{i,j} = [-1,1];
        else
            State_change{i,j} = [-1,-1];
        end
        S = State_change{i,j};
        State_change1(i,j) = S(1,1);
        State_change0(i,j) = S(1,2);        
    end
end
State_change1 = State_change1(:,1:2);
State_change0 = State_change0(:,1:2);
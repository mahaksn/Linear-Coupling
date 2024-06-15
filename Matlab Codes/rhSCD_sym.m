function [rhTable,stable]=rhSCD_sym(coeffVector)
coeffLength = length(coeffVector);
rhTableColumn = round(coeffLength/2);

rhTable(1,:) = coeffVector(1,1:2:coeffLength);

if (rem(coeffLength,2) ~= 0)
    rhTable(2,1:rhTableColumn - 1) = coeffVector(1,2:2:coeffLength);
else
    rhTable(2,:) = coeffVector(1,2:2:coeffLength);
end

%% Calculate Routh-Hurwitz table's rows
%  Set epss as a small value
epss = 0.01;
%  Calculate other elements of the table
for i = 3:coeffLength
   
    %  special case: row of all zeros
    if rhTable(i-1,:) == 0
        order = (coeffLength - i);
        cnt1 = 0;
        cnt2 = 1;
        for j = 1:rhTableColumn - 1
            rhTable(i-1,j) = (order - cnt1) * rhTable(i-2,cnt2);
            cnt2 = cnt2 + 1;
            cnt1 = cnt1 + 2;
        end
    end
    
    for j = 1:rhTableColumn - 1
        %  first element of upper row
        firstElemUpperRow = rhTable(i-1,1);
        
        %  compute each element of the table
        rhTable(i,j) = ((rhTable(i-1,1) * rhTable(i-2,j+1)) - ....
            (rhTable(i-2,1) * rhTable(i-1,j+1))) / firstElemUpperRow;
    end
    
    
    %  special case: zero in the first column
    if rhTable(i,1) == 0
        rhTable(i,1) = epss;
    end
end
%%  Compute number of right hand side poles(unstable poles)
%   Initialize unstable poles with zero
unstablePoles = 0;
%   Check change in signs
for i = 1:coeffLength - 1
    if sign(rhTable(i,1)) * sign(rhTable(i+1,1)) == -1
        unstablePoles = unstablePoles + 1;
    end
end
% fprintf('\n Routh-Hurwitz Table:\n')
rhTable;
%   Print the stability result on screen
if unstablePoles == 0
    % fprintf('~~~~~> it is a stable system! <~~~~~\n')
    stable=1;
else
    % fprintf('~~~~~> it is an unstable system! <~~~~~\n')
    stable=0;
end
% fprintf('\n Number of right hand side poles =%2.0f\n',unstablePoles)
end
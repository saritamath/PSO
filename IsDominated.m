function dominated = IsDominated(risks, returns)
    n = length(risks);
    dominated = false(n, 1);
    for i = 1:n
        for j = 1:n
            if i ~= j && risks(j) <= risks(i) && returns(j) >= returns(i)
                if risks(j) < risks(i) || returns(j) > returns(i)
                    dominated(i) = true;
                    break;
                end
            end
        end
    end
end

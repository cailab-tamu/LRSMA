function [intervals] = makeIntervals(data, threshold, gapMax)
data = find(data > threshold);
gap = find(diff(data) > gapMax);
intervals = [];
if length(gap) == 0 | length(unique(data)) < 3
    intervals = [min(data) max(data)];
else
    for i = 1:length(gap)
        if i == 1
            intervals = [intervals ; min(data) data(gap(i))];
            if i == length(gap)
                intervals = [intervals ; data(gap(i)+1) max(data)];
            end
        else
            if i == length(gap)
                intervals = [intervals ; data(gap(i)+1) max(data)];
            else
                intervals = [intervals; data(gap(i)+1) data(gap(i+1))];
            end
        end
    end
end
end
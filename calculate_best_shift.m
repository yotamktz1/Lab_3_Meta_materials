function [best_shift] = calculate_best_shift(set_1, set_2)

% this function must be used carefully, it's not idiot proof AT ALL!
% one should check first which set of points should be shifted and choose
% the resolution wisely

% max_shift = input("value for max shift");
resolution = input("resolution"); % the range of shifts, by number of frequency units

res_mat = zeros(2, resolution);

for i = 1:resolution
    diff_vec = set_1(1: end-i) - set_2(1+i :end);
    diff = mean(abs(diff_vec));
    res_mat(:,i) = [diff i];
end

[val,best_shift] = min( res_mat(1,:) ) ;
disp("best shift is : " +  best_shift + "it's value is" + val);

end

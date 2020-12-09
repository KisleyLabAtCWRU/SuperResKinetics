function points = change_point_detection(eff)
%% The main function to detect all the change points in a trace
sd = w1_noise(diff(eff))/sqrt(2);% estimate the noise level
points = recursion1(eff,sd,[], 0);% recursively detect all the change points
points = sort(points);
end

%% 
function points = recursion1(eff, sd, points, counter)
tau95 = 1;% the threshold to consider a spot as a change point
N = numel(eff);
if N < 2% only one point left in the segment, stop searching for the change point
    return
else
    llr = change_point_wavelet(eff,sd);
    [Z, k] = max(abs(llr));
    if Z > tau95
        counter = 0;
        points(end+1) = k;
        points1 = recursion1(eff(1:k), sd, [], counter);
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    elseif counter < 3% the parameter to dig in and find more short-lived transitions
        counter = counter +1;
        k = floor(numel(eff)/2);
        points1 = recursion1(eff(1:k), sd, [], counter);
        points2 = recursion1(eff(k+1:end), sd, [], counter);
        points = [points, points1, points2+k];
    else
        counter = 0;
        return
    end
end
end

%% combine the idea of Haar wavelet and change point method
function llr = change_point_wavelet(eff, sd)
N = numel(eff);
llr = zeros(size(eff));
for i = 1 : N-1
    I1 = mean(eff(1:i));
    I2 = mean(eff(i+1:end));
    llr(i) = (I2 - I1)/3/sd/sqrt(1/i+1/(N-i));
end
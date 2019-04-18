function main
close all;
load data;
load test_data;

X = 0:(length(test_data)-1);
hold on;
plot(data(:,1), data(:,2), '.-r');
plot(X, test_data, 'o-b');


end

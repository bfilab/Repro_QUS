clear; close all; clc;

data1 = normrnd(0,1,[10000,1]);
data2 = normrnd(0,1,[10000,1]);

subplot(2,1,1)
hist(data1,100);
xdim = xlim;
subplot(2,1,2)
hist(data2,100);
xlim([xdim])

[gCNR] = computeContrastMetrics(data1,data2);

title(gCNR)
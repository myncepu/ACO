clc; clear all; close all;
oF = @(x, y) -(x.^4 + 3*y.^4 - 0.2*cos(3*pi*x) - 0.4*cos(4*pi*y) + 0.6);
simpleACO(oF, 1, 100, 100, 0.9, 1, 0.05, [1, 1], [-1, -1])

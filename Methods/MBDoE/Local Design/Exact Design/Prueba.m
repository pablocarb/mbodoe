clc, clear all, close all

% Example 6.2
% 
% s = Sensivity(@(x, p)Exercise6_2(x, p), 0, [1.5, 3])
% 
% F = Fisher(@(x, p)Exercise6_2(x, p), [0, 1/3]' ,[1.5, 3], 1)
% 
% det(F)
% 
% [D_opt, j_D] = Fedorov(@(x,p)Exercise6_3(x, p), [1, 0.5, 1, 0.1], 1, transpose([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1, 1.2, 1.5, 10]), [], [], [], [], 0, Inf, [])

% Example 6.3

[D_opt, j_D_opt] = Detmax(@(x,p)Exercise6_3(x, p), [1, 1, 1, 0.1], 1, transpose([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.9, 1, 1.2, 1.5, 10]), 3, [], [], [], [], 0, Inf, [])

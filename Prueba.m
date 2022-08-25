clc; clear all; close all;



[D_opt, j_D_opt] = Detmax_disc_ga(@(x,p)Ex_model(x, p), [1, 1, 1, 0.1], 201, 12, 1, 1e-12, 6, [], [], [], [], 0, 20, []);
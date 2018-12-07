function [Rows,Cols] = optiSubPlotNum(AtomImgNum)
% Function to determine optimal number of rows and columns for subplot.
% Used to modularize code more effectively
%
% INPUTS:
%   AtomImgNum - Number of Atom images that will be plotted. This
%                determines how best to organize plots.
%
% OUTPUTS:
%   Rows - Number of subplot rows
%   Cols - Number of subplot rows
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if AtomImgNum == 1
    Rows = 1;
    Cols = 1;
elseif AtomImgNum == 2
    Rows = 1;
    Cols = 2;
elseif AtomImgNum <= 4
    Rows = 2;
    Cols = 2;
elseif AtomImgNum <= 9
    Rows = 3;
    Cols = 3;
elseif AtomImgNum <= 12
    Rows = 3;
    Cols = 4;
elseif AtomImgNum <= 16
    Rows = 4;
    Cols = 4;
elseif AtomImgNum <= 28
    Rows = 4;
    Cols = floor(AtomImgNum/4 - 0.01) + 1;
else
    Rows = 5;
    Cols = floor(AtomImgNum/5 - 0.01) + 1;
end
end
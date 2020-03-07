% A simplified little Matlab file to load in some data and 
% run the test for a monotonic relationship in expected returns on portfolios 
% formed by sorting on some variable, as described in:
%
% "Monotonicity in Asset Returns: New Tests with Applications to the Term
% Structure, the CAPM and Portfolio Sorts"
%
% Paper available at:
% http://www.economics.ox.ac.uk/members/andrew.patton/research.html
%
%
%  Andrew Patton 
%
% 29 April 2009

close all;clear all;clc

% THIS FILE was designed for Matlab beginners and I've tried to be very
% explicit about each step. If you already know a little bit of Matlab,
% then you can find more details in "sorts_example_code_1.m", which
% replicates all the results in the paper - you can then adapt that
% function to use on your own data.


% INSTRUCTIONS: read through the following, save your data as in Step 2,
% and enter the values requested in Step 3. Once you've reached the bottom
% of this file: save it, and then hit F5, which will run the code. 
% The MR p-value will appear in the Matlab "command window", and a Figure
% will be generated using your data.


% Please note: THIS FUNCTION ONLY WORKS FOR SORTS USING A SINGLE VARIABLE
% if you want to conduct a test based on sorts using 2 or more variables
% then please either see the more general example code
% (sorts_example_code_1.m) 


% STEP 1: Download the ZIP file "Patton_Timmermann_MR_test_toolbox.zip".
% Unzip the contents of this file into the directory on your computer where
% you have Matlab installed. For example: 
% "C:\Program Files\MATLAB\R2008a\toolbox\".
%
% This will give you all the Matlab functions you need to run this code.


% STEP 2: Get the returns from your portfolio sorts into a TxK matrix,
% where T is the number of time series observations and K is the number of portfolios. 
% Remove all text from this matrix (eg, column titles and row titles) so
% that you just have a matrix of numbers.
%
% Save this matrix as "data1.txt' in the "C:\" directory of your computer.


% STEP 3: answer the following questions:

% How many bootstrap replications would you like to use? (I
% recommend at least 1000)
bootreps = 1000;

% Do you want to test for an *increasing" relationship
% (direction=+1) or a *decreasing* relationship (direction=-1)?
direction = +1;

% What "block" size do you want to use? (this is related to how
% much serial correlation you think there is in your data). As a rough
% guide, for daily/monthy/quarterly/annual returns data I would recommend
% 10/6/3/2 as the block length, though you may wish to see how sensitive
% your results are to this choice.
block_length = 6;


rand_state = 1234;  % this just fixes the "seed" for the random number generator, to make sure you get the same answer every time you run this code.


load 'c:\data1.txt' -ascii;
[T,K] = size(data1);

MR_pvalue = MR_test_22_simple(data1,bootreps,direction,block_length,rand_state)

figure(1),plot((1:K)',mean(data1),'bo-'),...
    xlabel('Portfolio number'),ylabel('Average return'),...
    title(['Average returns on sorted portfolio returns, MR p-value = ',num2str(MR_pvalue,'%10.3f')]);


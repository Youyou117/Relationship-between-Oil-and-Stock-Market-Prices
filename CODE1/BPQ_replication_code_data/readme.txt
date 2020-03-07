Readme file for BPQ_2016_code_data.zip

For the paper:

Bollerslev, T., A. J. Patton, and R. Quaedvlieg, 2016, 
Exploiting the Errors: A Simple Approach for Improved Volatility Forecasting, 
Journal of Econometrics, 192, 1-18.


This ZIP file contains the following:


1) Nine Matlab files, containing code to replicate the tables in the above paper:

	a) BPQ2016_Replication_SP500.m, replicates all numbers relating to the S&P 500 index. (Takes around 3 minutes to run.) Need to adjust
	   lines 35 and 37 according to where you save the data files.

	b) BPQ2016_Replication_Stocks.m, replicates all numbers relating to the 27 individual stocks we consider. (Takes around 90 minutes to run.) Need
	   to adjust line 38 according to where you save the data files.

	c) The remaining Matlab files are helper functions that are called in the main two files. These are mostly from LeSage's Matlab Econometrics toolboxes.
	   http://www.spatial-econometrics.com/


2) Another ZIP file containing data for this paper:

	a) The folder "Data_5min" contains various daily measures of volatility (and quarticity) that used in the paper, based on 5-minute returns.

	b) The folder "Data_1min" contains various daily measures of volatility (and quarticity) that used in the paper, based on 1-minute returns.

	c) The folder "Market" contains two Excel spreadsheets containing a summary of the data, just for the S&P500 index.


3) A PDF file containing the tables that this code produces.





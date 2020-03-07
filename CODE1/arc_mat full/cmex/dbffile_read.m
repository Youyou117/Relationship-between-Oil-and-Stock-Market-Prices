function [data,vnames] = dbffile_read(filename)
% PURPOSE: reads an arcview (or other) dbffile and returns a structure containing
%          a data matrix and variable names vector
% -----------------------------------------------------
% USAGE: [data,vnames] = dbffile_read(filename)
% where: filename = an arcview (or dbf) file name without the extension
% ----------------------------------------------------------------------------------------
% Returns: data = an nobs x nvars matrix
%          vnames = an nvars x 1 vector of variable names
% ----------------------------------------------------------------------------------------
% NOTES: simply calls a c-mex function dbf_read
% compile with: mex dbf_read.c shapelib.c
% ----------------------------------------------------------------------------------------

% written by: James P. LeSage 12/2003
% University of Toledo
% Department of Economics
% Toledo, OH 43606
% jlesage@spatial-econometrics.com

[data,vnames] = dbf_read(filename);

vnames = strvcat(vnames);

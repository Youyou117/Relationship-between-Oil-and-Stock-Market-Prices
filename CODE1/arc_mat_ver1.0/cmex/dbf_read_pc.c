
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <malloc.h>
#include <mex.h>
#include "shapelib.h"


// *****************************************************
// mex dbf_read.c file
// compile with: mex dbf_read.c shapelib.c
// shapelib.c uses shapelib.h file as well
// this file is called by dbffile_read.m, a matlab file


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  char * filename;
  double *data;
  int   buflen,status;
  DBFHandle	hDBF;
  int	i, iRecord;
  char	szFormat[32];
  int		nWidth, nDecimals;
  int		bHeader = 0;
  int		bRaw = 0;
  int		bMultiLine = 0;
  char	szTitle[30];
  int nrecords, nfields, line;


  /* Check for proper number of arguments. */
    /* Check for proper number of arguments. */
  if (nrhs != 1) 
    mexErrMsgTxt("dbf_read: one filename input required.");
  else if (nlhs != 2) 
    mexErrMsgTxt("dbf_read: wrong # of output arguments.");

  /* Input must be a string. */
  if (mxIsChar(prhs[0]) != 1)
    mexErrMsgTxt("dbf_read: Input must be a string.");

  /* Input must be a row vector. */
  if (mxGetM(prhs[0]) != 1)
    mexErrMsgTxt("dbf_read: Input must be a row vector.");
    
  /* Get the length of the input string. */
  buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;

  /* Allocate memory for input string. */
  filename = mxCalloc(buflen, sizeof(char));

  /* Copy the string data from prhs[0] into a C string 
   * filename. If the string array contains several rows, 
   * they are copied, one column at a time, into one long 
   * string array. */
  status = mxGetString(prhs[0], filename, buflen);
  if (status != 0) 
    mexWarnMsgTxt("dbf_read: Not enough space. filename is truncated.");
    
/* -------------------------------------------------------------------- */
/*      Open the file.                                                  */
/* -------------------------------------------------------------------- */
    hDBF = DBFOpen( filename, "rb" );
    if( hDBF == NULL )
    {
    mexErrMsgTxt("dbf_read: Unable to open dbf file");
    }
    
/* -------------------------------------------------------------------- */
/*	If there is no data in this file let the user know.		*/
/* -------------------------------------------------------------------- */
    if( DBFGetFieldCount(hDBF) == 0 )
    {
    mexErrMsgTxt("dbf_read: There are no fields in this table!");
    }
    

 /* find out how much space to allocate for return entities */
 
 nrecords = DBFGetRecordCount(hDBF);
 nfields = DBFGetFieldCount(hDBF);

    /* Create matrices for the return arguments */
    plhs[0] = mxCreateDoubleMatrix(nrecords,nfields, mxREAL); // data matrix
    data = mxGetPr(plhs[0]);
    
 
/* -------------------------------------------------------------------- */
/*	Read all the records and dump data matrix  						*/
/* -------------------------------------------------------------------- */
    for( iRecord = 0; iRecord < nrecords; iRecord++ )
    {
	for( i = 0; i < nfields; i++ )
	{
    *(data + iRecord + i*nrecords) = DBFReadDoubleAttribute( hDBF, iRecord, i );                
	} // end of for loop over fields within a record
	
    } // end of for loop over records
    
    
  plhs[1] = mxCreateCellMatrix(nfields, 1);
  for(i=0; i<nfields; i++) {
     DBFFieldType	eType;
     eType = DBFGetFieldInfo( hDBF, i, szTitle, &nWidth, &nDecimals );
     mxSetCell(plhs[1], i, mxCreateString(szTitle));
  }
   
    DBFClose( hDBF );

}



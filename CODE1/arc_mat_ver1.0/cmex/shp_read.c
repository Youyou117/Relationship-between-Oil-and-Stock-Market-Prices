
#include <stdio.h>
#include <math.h>
#include <string.h>
//#include <malloc.h>
#include <mex.h>
#include "shapelib.h"


// *****************************************************
// mex shp_read.c file
// compile with: mex shp_read.c shapelib.c
// shapelib.c uses shapelib.h file as well
// this file is called by shpfile_read.m, a matlab file


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
  char * filename;
  int   buflen,status;
  double *cartex, *cartey, *xmin, *xmax, *ymin, *ymax, *nvertices, *nparts;
  int  n, np, cnt;
  SHPHandle   hSHP;
  int     nShapeType, nEntities, i, iPart; 
  int bValidate = 0,nInvalidCount=0;
  const char  *pszPlus;
  double  adfMinBound[4], adfMaxBound[4];
  double nan;
  double *xcc, *ycc, xc, yc, xct, yct;
  int     ring, ringPrev, ring_nVertices, rStart;
  double  Area, ringArea, Areat;


  /* Check for proper number of arguments. */
    /* Check for proper number of arguments. */
  if (nrhs != 1) 
    mexErrMsgTxt("shp_read: one filename input required.");
  else if (nlhs != 10) 
    mexErrMsgTxt("shp_read: wrong # of output arguments.");

  /* Input must be a string. */
  if (mxIsChar(prhs[0]) != 1)
    mexErrMsgTxt("shp_read: Input must be a string.");

  /* Input must be a row vector. */
  if (mxGetM(prhs[0]) != 1)
    mexErrMsgTxt("shp_read: Input must be a row vector.");
    
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
    mexWarnMsgTxt("shp_read: Not enough space. filename is truncated.");

    nan=mxGetNaN();

    hSHP = SHPOpen( filename, "rb" );

    if( hSHP == NULL )
    {
    mexErrMsgTxt("shp_read: Unable to open shape file");
    }

 /* find out how much space to allocate for return entities */
    SHPGetInfo( hSHP, &nEntities, &nShapeType, adfMinBound, adfMaxBound );

   
    n = 0;
    np = 0;
for( i = 0; i < nEntities; i++ )
    {
    int     j;
    SHPObject   *psShape;
    psShape = SHPReadObject( hSHP, i );
    n = n + psShape -> nVertices;
    if ( psShape -> nParts > 1) // add space for NaN separating each part
    n = n + psShape -> nParts;
    
  SHPDestroyObject( psShape );

} // end of for i loop

// add 1 space for each polygon to have an NaN terminator
   n = n + nEntities;
     
    /* Create matrices for the return arguments */
    plhs[0] = mxCreateDoubleMatrix(n,1, mxREAL); // cartex
    plhs[1] = mxCreateDoubleMatrix(n,1, mxREAL); // cartey
    plhs[2] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // xmin shape bounds
    plhs[3] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // xmax
    plhs[4] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // ymin
    plhs[5] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // ymax
    plhs[6] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // nvertices
    plhs[7] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // nparts
    plhs[8] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // x-centroid
    plhs[9] = mxCreateDoubleMatrix(nEntities,1, mxREAL); // y-centroid
    


    cartex = mxGetPr(plhs[0]);
    cartey = mxGetPr(plhs[1]);
    xmin = mxGetPr(plhs[2]);
    xmax = mxGetPr(plhs[3]);
    ymin = mxGetPr(plhs[4]);
    ymax = mxGetPr(plhs[5]);
    nvertices = mxGetPr(plhs[6]);
    nparts = mxGetPr(plhs[7]);
    xcc = mxGetPr(plhs[8]);
    ycc = mxGetPr(plhs[9]);

/* -------------------------------------------------------------------- */
/*  Skim over the list of shapes, printing all the vertices.    */
/* -------------------------------------------------------------------- */
    cnt = 0;
    for( i = 0; i < nEntities; i++ )
    {
    int     j;
        SHPObject   *psShape;

    psShape = SHPReadObject( hSHP, i );
    
    *(nvertices + i) = (double) psShape->nVertices;
    *(nparts + i) = (double) psShape->nParts;
    
     *(xmin + i) = (double) psShape->dfXMin;
     *(ymin + i) = (double) psShape->dfYMin;
     *(xmax + i) = (double) psShape->dfXMax;
     *(ymax + i) = (double) psShape->dfYMax;
     
   Areat = 0.0;
   xct = 0.0;
   yct = 0.0;
   
   /* for each ring in compound / complex object calc the ring cntrd        */
   
   ringPrev = psShape->nVertices;
   for ( ring = (psShape->nParts - 1); ring >= 0; ring-- ) {
     rStart = psShape->panPartStart[ring];
     ring_nVertices = ringPrev - rStart;
     
    xc = 0.0;
    yc = 0.0;
    ringArea = 0.0;
    RingCentroid_2d ( ring_nVertices, (double*) &(psShape->padfX [rStart]), 
    (double*) &(psShape->padfY [rStart]), &xc, &yc, &ringArea);  

     /* use Superposition of these rings to build a composite Centroid      */
     /* sum the ring centrds * ringAreas,  at the end divide by total area  */
     xct +=  xc * ringArea;
     yct +=  yc * ringArea; 
     Areat += ringArea; 
     ringPrev = rStart;
    }    

     /* hold on the division by AREA until were at the end                  */
     xc = xct / Areat;
     yc = yct / Areat;
    
   *(xcc + i) = xc;
   *(ycc + i) = yc;
    
        for( j = 0, iPart = 1; j < psShape->nVertices; j++ )
    {
        const char  *pszPartType = "";

        if( j == 0 && psShape->nParts > 0 )
                pszPartType = SHPPartTypeName( psShape->panPartType[0] );
            
        if( iPart < psShape->nParts && psShape->panPartStart[iPart] == j )
        {
                pszPartType = SHPPartTypeName( psShape->panPartType[iPart] );
        iPart++;
        // starting a new part, so put an NaN separator here
       *(cartex + cnt) = nan;
       *(cartey + cnt) = nan;
        cnt = cnt + 1;
       *(cartex + cnt) = psShape->padfX[j];
       *(cartey + cnt) = psShape->padfY[j];
        cnt = cnt + 1;
        }
        else {
      *(cartex + cnt) = psShape->padfX[j];
      *(cartey + cnt) = psShape->padfY[j];
      cnt = cnt + 1;
      }

    } /* end of for j loop */
    
      *(cartex + cnt) = nan;
      *(cartey + cnt) = nan;
      cnt = cnt + 1;
    SHPDestroyObject( psShape );
    } // end of for i loop
    
    
    SHPClose( hSHP );

}

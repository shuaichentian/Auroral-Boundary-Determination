
// Author: Xavier Bresson (xbresson at math.ucla.edu)
// Last version: Oct 06 2008
// Name: sdf_mex.c
// Description: This function computes the signed distance function of the zero level set
// of the input function "pfDistSDF". The re-distancing the level set function is done with
// the Fast Marching Method of Adalsteinsson and Sethian [D. Adalsteinsson and J. Sethian, "A Fast Level
// Set Method for Propagating Interfaces," Journal of Computational Physics, vol. 118, pp. 269-277,1995.]

// Compilation (run in matlab):
// mex -v -g sdf_mex.c 
// mex sdf_mex.c 



#include <mex.h>
#include <math.h>
#include <time.h>

#define TRUE 1

#define YES 0
#define NO 1

#define MAX(a,b) ( a > b ? a : b )
#define MIN(a,b) ( a < b ? a : b )

#define MAX5(a,b,c,d,e) MAX( MAX(MAX(a,b), MAX(c,d)), e )
#define MIN5(a,b,c,d,e) MIN( MIN(MIN(a,b), MIN(c,d)), e )

#define SIGN(x) ( x > 0.0 ? 1.0 : -1.0 )

#define SWAP(a,b,tmp) tmp=a; a=b; b=tmp

#define ABS(x) ( x >= 0.0 ? x : -x )
#define SQR(x) (x)*(x)

#define INF 1024.0
#define EPS 0.001




/****************************************/
/**  Declaration of structures         **/
/****************************************/

typedef struct{
  float *pfPhiValue;
  int   *pix;
  int   *piy;
  int   *piFlagNarrowBand;
  int   iLength;
} sNarrowBand;


/****************************************/



/****************************************/
/**  Definition of static procedures   **/
/****************************************/


/****************************************/
void vAddElement2NarrowBand
(
 sNarrowBand  *psNarBand,
 float        fPhiValue,
 int          ix,
 int          iy
)

{
  
  int    current, parent, iTmp;
  float  fTmp;

  psNarBand->pfPhiValue[psNarBand->iLength] = fPhiValue;
  psNarBand->pix[psNarBand->iLength] = ix;
  psNarBand->piy[psNarBand->iLength] = iy;

  current = psNarBand->iLength;
  psNarBand->iLength++;
  
  parent = (current - 1)/2;
  while (current != 0)
  {
      if (ABS(psNarBand->pfPhiValue[current]) < ABS(psNarBand->pfPhiValue[parent]))
      {
          SWAP(psNarBand->pfPhiValue[current],psNarBand->pfPhiValue[parent],fTmp);
          SWAP(psNarBand->pix[current],psNarBand->pix[parent],iTmp);
          SWAP(psNarBand->piy[current],psNarBand->piy[parent],iTmp);
          current = parent;
          parent = (current-1)/2;
      }
      else
          break;
  }
  
}
/****************************************/


/****************************************/
void vRemoveFirstElementOfNarrowBand
(
 sNarrowBand  *psNarBand
)

{

  int    best, c1, c2, iTmp, current;
  float  fTmp;


  psNarBand->pfPhiValue[0] = psNarBand->pfPhiValue[psNarBand->iLength-1];
  psNarBand->pix[0] = psNarBand->pix[psNarBand->iLength-1];
  psNarBand->piy[0] = psNarBand->piy[psNarBand->iLength-1];
  psNarBand->iLength--;
 
  current = 0;
  c1 = 1;
  c2 = 2;
  while (c1 < psNarBand->iLength)
  {
      if (c2 >= psNarBand->iLength)
          best = c1;
      else
      {
          if (ABS(psNarBand->pfPhiValue[c1]) < ABS(psNarBand->pfPhiValue[c2]))
              best = c1;
          else
              best = c2;
      }
      if (ABS(psNarBand->pfPhiValue[best]) < ABS(psNarBand->pfPhiValue[current]))
      {
          SWAP(psNarBand->pfPhiValue[best],psNarBand->pfPhiValue[current],fTmp);
          SWAP(psNarBand->pix[best],psNarBand->pix[current],iTmp);
          SWAP(psNarBand->piy[best],psNarBand->piy[current],iTmp);
          current = best;
          c1 = 2*current + 1;
          c2 = c1 + 1;
      }
      else
          break;
  }
  
}
/****************************************/


/****************************************/
void vUpdatePixelFromNarrowBand
(
 sNarrowBand  *psNarBand,
 float        fPhiValue,
 int          ix,
 int          iy
)

{

  int    current, parent, iTmp, i;
  float  fTmp;


  for(i=0; i< psNarBand->iLength; i++)
      if ( (psNarBand->pix[i] == ix) && (psNarBand->piy[i] == iy) )
          current = i;
  
  psNarBand->pfPhiValue[current] = fPhiValue;
  
  parent = (current - 1)/2;
  while (current != 0)
  {
      if (ABS(psNarBand->pfPhiValue[current]) < ABS(psNarBand->pfPhiValue[parent]))
      {
          SWAP(psNarBand->pfPhiValue[current],psNarBand->pfPhiValue[parent],fTmp);
          SWAP(psNarBand->pix[current],psNarBand->pix[parent],iTmp);
          SWAP(psNarBand->piy[current],psNarBand->piy[parent],iTmp);
          current = parent;
          parent = (current-1)/2;
      }
      else
          break;
  }
  
}
/****************************************/


/****************************************/
int iFM2D_SDF
(
 float  *pfDistSDF,
 float  *pfNewSDF,
 int    iDim[2],
 float  fMaxDistanceComputations,
 int    iNbIterZLS
 )

{

  float   fPhi00, fPhi0p, fPhi0m, fPhip0, fPhim0, fPhi1, fPhi2;
  float   fx1, fThreshold;
  float   fDeltaT, fDxp, fDxm, fDyp, fDym, fGrad, fs, fPhiValue;
  float   *pfDistSDFt, *pfNewSDFt;
  int     iNy, iNx, iy, ix, iy2, ix2, iCpt, iCptNbIter, iLengthNarrowBandTemp;
  int     *piBurntPixels, iNy2, iNx2;
  sNarrowBand   sNarBand, sNarBandTemp;


  /* Get the size of the image */
  iNy = iDim[0]; 
  iNx = iDim[1]; 

  iNy2 = iNy + 2;
  iNx2 = iNx + 2;



  /*****************************************************************/
  /* Memory allocations */

  pfDistSDFt = (float *) calloc( iNy2* iNx2,sizeof(float) );
  if (!pfDistSDFt)
    {
      mexPrintf("Memory allocation failure for pfDistSDFt\n");
      return(0);
    }

  pfNewSDFt = (float *) calloc( iNy2* iNx2,sizeof(float) );
  if (!pfNewSDFt)
    {
      mexPrintf("Memory allocation failure for pfNewSDFt\n");
      return(0);
    }

  sNarBand.pfPhiValue = (float *) calloc( iNy2* iNx2,sizeof(float) );
  if (!sNarBand.pfPhiValue)
    {
      mexPrintf("Memory allocation failure for sNarBand.pfPhiValue\n");
      return(0);
    }

  sNarBand.pix = (int *) calloc( iNy2* iNx2,sizeof(int) );
  if (!sNarBand.pix)
    {
      mexPrintf("Memory allocation failure for sNarBand.pix\n");
      return(0);
    }

  sNarBand.piy = (int *) calloc( iNy2* iNx2,sizeof(int) );
  if (!sNarBand.piy)
    {
      mexPrintf("Memory allocation failure for sNarBand.piy\n");
      return(0);
    }

  sNarBand.piFlagNarrowBand = (int *) calloc( iNy2* iNx2,sizeof(int) );
  if (!sNarBand.piFlagNarrowBand)
    {
      mexPrintf("Memory allocation failure for sNarBand.piFlagNarrowBand\n");
      return(0);
    }

  sNarBand.iLength = 0;

  piBurntPixels = (int *) calloc( iNy2* iNx2,sizeof(int) );
  if (!piBurntPixels)
    {
      mexPrintf("Memory allocation failure for piBurntPixels\n");
      return(0);
    }


  /*****************************************************************/

  
  for(ix=1; ix< iNx2-1; ix++)
      for(iy=1; iy< iNy2-1; iy++)
  {
      piBurntPixels[ix*iNy2+ iy] = 0;
      sNarBand.piFlagNarrowBand[ix*iNy2+ iy] = 0;
      
      if ( pfDistSDF[(ix-1)*iNy + iy-1] >= 0.0 )
          pfNewSDFt[ix*iNy2 + iy] = INF;
      else
          pfNewSDFt[ix*iNy2 + iy] = -INF;
      
      pfDistSDFt[ix*iNy2 + iy] = pfDistSDF[(ix-1)*iNy + iy-1];
      }
  
  
  for(ix=0; ix< iNx2; ix++)
  {
      piBurntPixels[ix*iNy2] = 1;
      sNarBand.piFlagNarrowBand[ix*iNy2] = 1;
      
      piBurntPixels[ix*iNy2+ iNy2-1] = 1;
      sNarBand.piFlagNarrowBand[ix*iNy2+ iNy2-1] = 1;
      
      if ( pfDistSDFt[ix*iNy2+ 1]>= 0.0 )
      {
          pfDistSDFt[ix*iNy2] = INF;
          pfNewSDFt[ix*iNy2] = INF;
      }
      else
      {
          pfDistSDFt[ix*iNy2] = -INF;
          pfNewSDFt[ix*iNy2] = -INF;
      }
      if ( pfDistSDFt[ix*iNy2+ iNy2-2]>= 0.0 )
      {
          pfDistSDFt[ix*iNy2+ iNy2-1] = INF;
          pfNewSDFt[ix*iNy2+ iNy2-1] = INF;
      }
      else
      {
          pfDistSDFt[ix*iNy2+ iNy2-1] = -INF;
          pfNewSDFt[ix*iNy2+ iNy2-1] = -INF;
      }
  }
  
  
  for(iy=0; iy< iNy2; iy++)
  {
      piBurntPixels[iy] = 1;
      sNarBand.piFlagNarrowBand[iy] = 1;
      
      piBurntPixels[(iNx2-1)*iNy2+ iy] = 1;
      sNarBand.piFlagNarrowBand[(iNx2-1)*iNy2+ iy] = 1;
      
      if ( pfDistSDFt[iNy2+ iy]>= 0.0 )
      {
          pfDistSDFt[iy] = INF;
          pfNewSDFt[iy] = INF;
      }
      else
      {
          pfDistSDFt[iy] = -INF;
          pfNewSDFt[iy] = -INF;
      }
      if ( pfDistSDFt[(iNx2-2)*iNy2+ iy]>= 0.0 )
      {
          pfDistSDFt[(iNx2-1)*iNy2+ iy] = INF;
          pfNewSDFt[(iNx2-1)*iNy2+ iy] = INF;
      }
      else
      {
          pfDistSDFt[(iNx2-1)*iNy2+ iy] = -INF;
          pfNewSDFt[(iNx2-1)*iNy2+ iy] = -INF;
      }
  }
  



  fDeltaT = 0.6;
  for(ix=1; ix< iNx2-1; ix++)
      for(iy=1; iy< iNy2-1; iy++)
  {
      
      fPhi00 = pfDistSDFt[ix*iNy2 + iy];
      fPhip0 = pfDistSDFt[(ix+1)*iNy2 + iy];
      fPhim0 = pfDistSDFt[(ix-1)*iNy2 + iy];
      fPhi0p = pfDistSDFt[ix*iNy2 + iy+1];
      fPhi0m = pfDistSDFt[ix*iNy2 + iy-1];
      
      if ( MAX5(fPhi00,fPhip0,fPhim0,fPhi0p,fPhi0m)>= 0.0  &&  MIN5(fPhi00,fPhip0,fPhim0,fPhi0p,fPhi0m)< 0.0 )
      {
          
          fDxp = pfDistSDFt[(ix+1)*iNy2 + iy] - pfDistSDFt[ix*iNy2 + iy];
          fDxm = pfDistSDFt[ix*iNy2 + iy] -     pfDistSDFt[(ix-1)*iNy2 + iy];
          fDyp = pfDistSDFt[ix*iNy2 + iy+1] -   pfDistSDFt[ix*iNy2 + iy];
          fDym = pfDistSDFt[ix*iNy2 + iy] -     pfDistSDFt[ix*iNy2 + iy-1];
          if (fPhi00 >= 0.0)
              fGrad = sqrt( SQR( MAX(-MIN(fDxp,0.0), MAX(fDxm,0.0)) ) + SQR( MAX(-MIN(fDyp,0.0), MAX(fDym,0.0)) ) );
          else
              fGrad = sqrt( SQR( MAX(MAX(fDxp,0.0), -MIN(fDxm,0.0)) ) + SQR( MAX(MAX(fDyp,0.0), -MIN(fDym,0.0)) ) );
          fs = fPhi00/ sqrt(SQR(fPhi00) + SQR(fGrad) );
          
          pfNewSDFt[ix*iNy2+ iy] = pfDistSDFt[ix*iNy2+ iy] - fDeltaT* fs* (fGrad - 1.0);
          
          
          vAddElement2NarrowBand(&sNarBand,pfNewSDFt[ix*iNy2+ iy],ix,iy);
          sNarBand.piFlagNarrowBand[ix*iNy2+ iy] = 1;
          piBurntPixels[ix*iNy2+ iy] = 1;
          
      }
      
      }


  iLengthNarrowBandTemp = sNarBand.iLength;


  sNarBandTemp.pfPhiValue = (float *) calloc( iLengthNarrowBandTemp ,sizeof(float) );
  if (!sNarBandTemp.pfPhiValue)
    {
      mexPrintf("Memory allocation failure for sNarBandTemp.pfPhiValue\n");
      return(0);
    }

  sNarBandTemp.pix = (int *) calloc( iLengthNarrowBandTemp ,sizeof(int) );
  if (!sNarBandTemp.pix)
    {
      mexPrintf("Memory allocation failure for sNarBandTemp.pix\n");
      return(0);
    }

  sNarBandTemp.piy = (int *) calloc( iLengthNarrowBandTemp ,sizeof(int) );
  if (!sNarBandTemp.piy)
    {
      mexPrintf("Memory allocation failure for sNarBandTemp.piy\n");
      return(0);
    }

  sNarBandTemp.iLength = 0;

 
//   iNbIterZLS = 0;
//   mexPrintf("iNbIterZLS = %i\n",iNbIterZLS);

  for(iCpt=0; iCpt< iLengthNarrowBandTemp; iCpt++)
    {
      ix = sNarBand.pix[iCpt];
      iy = sNarBand.piy[iCpt];
      pfDistSDFt[ix*iNy2 + iy] = pfNewSDFt[ix*iNy2 + iy];
    }


  for(iCptNbIter=0; iCptNbIter< iNbIterZLS; iCptNbIter++)
    {
    
//         mexPrintf("Here\n");

        
      for(iCpt=0; iCpt< iLengthNarrowBandTemp; iCpt++)
	{

	  ix = sNarBand.pix[0];
	  iy = sNarBand.piy[0];
	  vRemoveFirstElementOfNarrowBand(&sNarBand);
	  fPhi00 = pfDistSDFt[ix*iNy2 + iy];
	  fDxp = pfDistSDFt[(ix+1)*iNy2 + iy] - pfDistSDFt[ix*iNy2 + iy];
	  fDxm = pfDistSDFt[ix*iNy2 + iy] -     pfDistSDFt[(ix-1)*iNy2 + iy];
	  fDyp = pfDistSDFt[ix*iNy2 + iy+1] -   pfDistSDFt[ix*iNy2 + iy];
	  fDym = pfDistSDFt[ix*iNy2 + iy] -     pfDistSDFt[ix*iNy2 + iy-1];
	  if (fPhi00 >= 0.0)
	    fGrad = sqrt( SQR( MAX(-MIN(fDxp,0.0), MAX(fDxm,0.0)) ) + SQR( MAX(-MIN(fDyp,0.0), MAX(fDym,0.0)) ) );
	  else
	    fGrad = sqrt( SQR( MAX(MAX(fDxp,0.0), -MIN(fDxm,0.0)) ) + SQR( MAX(MAX(fDyp,0.0), -MIN(fDym,0.0)) ) );
	  fs = fPhi00/ sqrt(SQR(fPhi00) + SQR(fGrad) );

	  pfNewSDFt[ix*iNy2 + iy] = pfDistSDFt[ix*iNy2 + iy] - fDeltaT* fs* (fGrad - 1.0);
	  
	  vAddElement2NarrowBand(&sNarBandTemp,pfNewSDFt[ix*iNy2+ iy],ix,iy);

	}

      for(iCpt=0; iCpt< iLengthNarrowBandTemp; iCpt++)
	{

	  ix = sNarBandTemp.pix[0];
	  iy = sNarBandTemp.piy[0];
	  fPhiValue = sNarBandTemp.pfPhiValue[0];
	  vRemoveFirstElementOfNarrowBand(&sNarBandTemp);
	  pfDistSDFt[ix*iNy2 + iy] = fPhiValue;
	  vAddElement2NarrowBand(&sNarBand,pfNewSDFt[ix*iNy2+ iy],ix,iy);

	}

    }


  fThreshold = 2.0;
  while ( sNarBand.iLength > 0 )
    {


      /* Take the point in "narrowband" that has the smallest value */
      ix = sNarBand.pix[0];
      iy = sNarBand.piy[0];
     
      /* Add the point to "burnt" and remove it from "narrowband" */
      piBurntPixels[ix*iNy2+ iy] = 1;
      vRemoveFirstElementOfNarrowBand(&sNarBand);
      sNarBand.piFlagNarrowBand[ix*iNy2+ iy] = 0;

      for(ix2=ix-1; ix2<= ix+1; ix2++)
	for(iy2=iy-1; iy2<= iy+1; iy2++)
	  if ( ix2>0  &&  ix2<iNx2-1  &&  iy2>0  &&  iy2<iNy2-1 )
	    if ( piBurntPixels[ix2*iNy2+ iy2] == 0 )
	      {

		if ( pfNewSDFt[ix2*iNy2+ iy2]> 0.0 )
		  {
		    fPhi1 = MIN(pfNewSDFt[(ix2+1)*iNy2+ iy2],pfNewSDFt[(ix2-1)*iNy2+ iy2]);
		    fPhi2 = MIN(pfNewSDFt[ix2*iNy2+ iy2+1],  pfNewSDFt[ix2*iNy2+ iy2-1]);
		  }
		else
		  {
		    fPhi1 = MIN(ABS(pfNewSDFt[(ix2+1)*iNy2+ iy2]),ABS(pfNewSDFt[(ix2-1)*iNy2+ iy2]));
		    fPhi2 = MIN(ABS(pfNewSDFt[ix2*iNy2+ iy2+1]),  ABS(pfNewSDFt[ix2*iNy2+ iy2-1]));
		  }

		if ( fPhi1==INF || fPhi2==INF )
		  {
		    if ( pfNewSDFt[ix2*iNy2+ iy2]> 0.0 )
		      pfNewSDFt[ix2*iNy2+ iy2] = MIN(fPhi1,fPhi2) + 1.0;
		    else
		      pfNewSDFt[ix2*iNy2+ iy2] = -MIN(fPhi1,fPhi2) - 1.0;
		  }
		else
		  {
		    fx1 = SQR(fPhi1- fPhi2);
		    if( fx1 < fThreshold )
		      if ( pfNewSDFt[ix2*iNy2+ iy2]> 0.0 )
			pfNewSDFt[ix2*iNy2+ iy2] = 0.5* ( fPhi1+ fPhi2+ sqrt(2.0- fx1) );
		      else
			pfNewSDFt[ix2*iNy2+ iy2] = -0.5* ( fPhi1+ fPhi2+ sqrt(2.0- fx1) );
		    else
		      if ( pfNewSDFt[ix2*iNy2+ iy2]> 0.0 )
			pfNewSDFt[ix2*iNy2+ iy2] = MIN(fPhi1,fPhi2) + 1.0;
		      else
			pfNewSDFt[ix2*iNy2+ iy2] = -MIN(fPhi1,fPhi2) - 1.0;
		  }



		if ( ABS(pfNewSDFt[ix2*iNy2+ iy2]) <= fMaxDistanceComputations )
		  {

		    if ( sNarBand.piFlagNarrowBand[ix2*iNy2+ iy2] == 0 )
		      {
			sNarBand.piFlagNarrowBand[ix2*iNy2+ iy2] = 1;
			vAddElement2NarrowBand(&sNarBand,pfNewSDFt[ix2*iNy2+ iy2],ix2,iy2);
		      }
		    else
		      vUpdatePixelFromNarrowBand(&sNarBand,pfNewSDFt[ix2*iNy2+ iy2],ix2,iy2);

		  }


	      }
    
    }



  for(ix=1; ix< iNx2-1; ix++)
    for(iy=1; iy< iNy2-1; iy++)
      if ( ABS(pfNewSDFt[ix*iNy2+ iy]) == INF )
	pfNewSDF[(ix-1)*iNy+ iy-1] = pfDistSDF[(ix-1)*iNy+ iy-1];
      else
	pfNewSDF[(ix-1)*iNy+ iy-1] = pfNewSDFt[ix*iNy2+ iy];
	   


  /*****************************************************************/
  /* Free memory  */

  free((float*) pfDistSDFt);
  free((float*) pfNewSDFt);
  free((float*) sNarBand.pfPhiValue);
  free((int*)   sNarBand.pix); 
  free((int*)   sNarBand.piy); 
  free((int*)   sNarBand.piFlagNarrowBand); 
  free((int*)   piBurntPixels);

  free((float*) sNarBandTemp.pfPhiValue);
  free((int*)   sNarBandTemp.pix);
  free((int*)   sNarBandTemp.piy);

  /*****************************************************************/

  return(1);
 
}
/****************************************/


/****************************************/




/*********************************************************************************************/
/*********************************************************************************************/
extern void mexFunction(int iNbOut, mxArray *pmxOut[],
		 int iNbIn, const mxArray *pmxIn[])
{

  /* iNbOut: number of outputs
     pmxOut: array of pointers to output arguments */

  /* iNbIn: number of inputs
     pmxIn: array of pointers to input arguments */

  
   float   *pfDistSDF, *pfNewSDF, *pfMaxDistanceComputations, fMaxDistanceComputations;
   int     iNbRows, iNbCols, iNx, iNy, i, iDim[3], iNdim;
   int     *piDisplay, iDisplay, *piNbIterZLS, iNbIterZLS;
   time_t  start_time, end_time;


  start_time = clock();


  /* Get the initial signed distance function (SDF) */
  pfDistSDF = mxGetData(pmxIn[0]);
  iNbRows = mxGetM(pmxIn[0]);     /* number of rows */
  iNbCols = mxGetN(pmxIn[0]);     /* number of columns */
  iNx = iNbCols;                  /* = Nx */
  iNy = iNbRows;                  /* = Ny */


  /* Get the maximal distance of computations */
  pfMaxDistanceComputations = mxGetData(pmxIn[1]);
  fMaxDistanceComputations = *pfMaxDistanceComputations;


  /* Get the  */
  piNbIterZLS = mxGetData(pmxIn[2]);
  iNbIterZLS = *piNbIterZLS;

  
  /* Get the  */
  piDisplay = mxGetData(pmxIn[3]);
  iDisplay = *piDisplay;

  


  iNdim = 2;
  iDim[0] = iNy;
  iDim[1] = iNx;
  
  pmxOut[0] = mxCreateNumericArray(iNdim,(const int*)iDim,mxSINGLE_CLASS,mxREAL);
  pfNewSDF = mxGetData(pmxOut[0]);


	
  if ( !iFM2D_SDF(pfDistSDF,pfNewSDF,iDim,fMaxDistanceComputations,iNbIterZLS) )
    {
      mexPrintf("\n\nError in SDF_FM: return an Array with null elements\n");
      for (i=0; i< iDim[0]* iDim[1]; i++)  pfNewSDF[i] = 0.0;
    }
  else
    {
      end_time = clock();      
	  if (iDisplay==YES) mexPrintf("\nCompute Signed Distance Function= %.5f sec\n\n",difftime(end_time,start_time)/1000);	
    }

}
/*********************************************************************************************/
/*********************************************************************************************/

/**************************************** End of file ****************************************/

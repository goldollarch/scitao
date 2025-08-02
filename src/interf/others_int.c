/*
scitao - A toolbox based on Scilab/Scicos environment for the simulation 
of wave optics, especially for the simulation of adaptive optics .

Copyright (c) 2005-2006 IAPCM, Beijing, China.  Written by
Chen jingyuan.  For comments or questions about this software,
please contact the author at jingyuan_chen@yahoo.com.cn.

This program is free software; you can redistribute it and/or modify it
under the terms of the GNU General Public License as  published by the
Free Software Foundation; either version 2 of the License, or (at your
option) any later version.

This program is provided "as is" and distributed in the hope that it
will be useful, but WITHOUT ANY WARRANTY; without even the implied
warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation,
Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA.
*/

#include "fitsio.h"
#include "int_optics.h"

int list_fits_head_int( char *fname )
{
	int m1, n1, l1, minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;
	
    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char card[FLEN_CARD];   /* Standard string lengths defined in fitsio.h */
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, nkeys, ii;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "c", &m1, &n1, &l1);

    if (!fits_open_file(&fptr, cstk(l1), READONLY, &status))
    {
      fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

      /* List only a single header if a specific extension was given */ 
      if (hdupos != 1 || strchr(cstk(l1), '[')) single = 1;

      for (; !status; hdupos++)  /* Main loop through each extension */
      {
        fits_get_hdrspace(fptr, &nkeys, NULL, &status); /* get # of keywords */

        sciprint("Header listing for HDU #%d:\n", hdupos);

        for (ii = 1; ii <= nkeys; ii++) { /* Read and print each keywords */

           if (fits_read_record(fptr, ii, card, &status))break;
           sciprint("%s\n", card);
        }
        sciprint("END\n\n");  /* terminate listing with END */

        if (single) break;  /* quit if only listing a single header */

        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try to move to next HDU */
      }

      if (status == END_OF_FILE)  status = 0; /* Reset after normal error */

      fits_close_file(fptr, &status);
    }

	if (status) sciprint("operate error !"); /* print any error message */

	return 0;
}

int list_fits_structure_int( char *fname )
{
	int m1, n1, l1, minlhs=1, maxlhs=1, minrhs=1, maxrhs=1;

    fitsfile *fptr;         /* FITS file pointer, defined in fitsio.h */
    char keyname[FLEN_KEYWORD], colname[FLEN_VALUE], coltype[FLEN_VALUE];
    int status = 0;   /* CFITSIO status value MUST be initialized to zero! */
    int single = 0, hdupos, hdutype, bitpix, naxis, ncols, ii;
    long naxes[10], nrows;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "c", &m1, &n1, &l1);

    if (!fits_open_file(&fptr, cstk(l1), READONLY, &status))
    {
      fits_get_hdu_num(fptr, &hdupos);  /* Get the current HDU position */

      /* List only a single structure if a specific extension was given */ 
      if (strchr(cstk(l1), '[') || strchr(cstk(l1), '+')) single++;

      for (; !status; hdupos++)   /* Main loop for each HDU */
      {
        fits_get_hdu_type(fptr, &hdutype, &status);  /* Get the HDU type */

        sciprint("\nHDU #%d  ", hdupos);
        if (hdutype == IMAGE_HDU)   /* primary array or image HDU */
        {
          fits_get_img_param(fptr, 10, &bitpix, &naxis, naxes, &status);

          sciprint("Array:  NAXIS = %d,  BITPIX = %d\n", naxis, bitpix);
          for (ii = 0; ii < naxis; ii++)
            sciprint("   NAXIS%d = %ld\n",ii+1, naxes[ii]);  
        }
        else  /* a table HDU */
        {
          fits_get_num_rows(fptr, &nrows, &status);
          fits_get_num_cols(fptr, &ncols, &status);

          if (hdutype == ASCII_TBL)
            sciprint("ASCII Table:  ");
          else
            sciprint("Binary Table:  ");

          sciprint("%d columns x %ld rows\n", ncols, nrows);
          sciprint(" COL NAME             FORMAT\n");

          for (ii = 1; ii <= ncols; ii++)
          {
            fits_make_keyn("TTYPE", ii, keyname, &status); /* make keyword */
            fits_read_key(fptr, TSTRING, keyname, colname, NULL, &status);
            fits_make_keyn("TFORM", ii, keyname, &status); /* make keyword */
            fits_read_key(fptr, TSTRING, keyname, coltype, NULL, &status);

            sciprint(" %3d %-16s %-16s\n", ii, colname, coltype);
          }
        }

        if (single) break;  /* quit if only listing a single HDU */

        fits_movrel_hdu(fptr, 1, NULL, &status);  /* try move to next ext */
      }

      if (status == END_OF_FILE) status = 0; /* Reset normal error */
      fits_close_file(fptr, &status);
    }

	if (status) sciprint("operate error !"); /* print any error message */

	return 0;
}

int get_fits_image_int( char *fname )
{
	int m1,  n1,  l1,  m2=1,  n2=2,  m3,  n3,
		minlhs=2, maxlhs=2, minrhs=1, maxrhs=1;

    fitsfile *fptr;  /* FITS file pointer */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int naxes[2], hdutype, naxis, *out_axes;
    long totpix, fpixel[2];
    double *pix;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "c", &m1, &n1, &l1);

    if ( !fits_open_file(&fptr, cstk(l1), READONLY, &status) )
    {
      if (fits_get_hdu_type(fptr, &hdutype, &status) || hdutype != IMAGE_HDU) { 
        sciprint("Error: this program only works on images, not tables\n");
        return(1);
      }

      fits_get_img_dim(fptr, &naxis, &status);
      fits_get_img_size(fptr, 2, naxes, &status);

      if (status || naxis != 2) { 
        sciprint("Error: NAXIS = %d.  Only 2-D images are supported.\n", naxis);
        return(1);
      }

      totpix = naxes[0] * naxes[1];
      pix = (double *) malloc( totpix * sizeof(double)); /* memory  */
      if (pix == NULL) {
        sciprint("Memory allocation error\n");
        return(1);
      }

      fpixel[0] = 1; fpixel[1] = 1;
	  fits_read_pix(fptr, TDOUBLE, fpixel, totpix, 0, pix, 0, &status);
     
      fits_close_file(fptr, &status);

    }

    if (status)  sciprint("operate error !"); /* print any error message */    

	out_axes =(int*)malloc(2*sizeof(int));
	out_axes[0]=naxes[0]; out_axes[1]=naxes[1];
	CreateVarFromPtr(2, "i", &m2, &n2, &out_axes);   
	LhsVar(1) = 2;

	m3=naxes[0]; n3=naxes[1];
	CreateVarFromPtr(3, "d", &m3, &n3, &pix);   
	LhsVar(2) = 3;

	return 0;
}

int fits_image_info_int( char *fname )
{
	int m1,  n1,  l1,  m2=1,  n2=2,  m3=1,  n3=4,
		minlhs=2, maxlhs=2, minrhs=1, maxrhs=1;

    fitsfile *fptr;  /* FITS file pointer */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int naxes[2], hdutype, naxis, ii, *out_axes=NULL;
    long totpix, fpixel[2];

	double *pix=NULL, *outptr=NULL;
    double sum = 0., meanval = 0., minval = 1.E33, maxval = -1.E33;

	CheckRhs(minrhs,maxrhs) ;
	CheckLhs(minlhs,maxlhs) ;

	GetRhsVar(1, "c", &m1, &n1, &l1);

    if ( !fits_open_file(&fptr, cstk(l1), READONLY, &status) )
    {
      if (fits_get_hdu_type(fptr, &hdutype, &status) || hdutype != IMAGE_HDU) { 
        sciprint("Error: this program only works on images, not tables\n");
        return(1);
      }

      fits_get_img_dim(fptr, &naxis, &status);
      fits_get_img_size(fptr, 2, naxes, &status);

      if (status || naxis != 2) { 
        sciprint("Error: NAXIS = %d.  Only 2-D images are supported.\n", naxis);
        return(1);
      }

      pix = (double *) malloc(naxes[0] * sizeof(double)); /* memory for 1 row */

      if (pix == NULL) {
        sciprint("Memory allocation error\n");
        return(1);
      }

      totpix = naxes[0] * naxes[1];

      fpixel[0] = 1;  /* read starting with first pixel in each row */

      /* process image one row at a time; increment row # in each loop */
      for (fpixel[1] = 1; fpixel[1] <= naxes[1]; fpixel[1]++)
      {  
         /* give starting pixel coordinate and number of pixels to read */
         if (fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0],0, pix,0, &status))
            break;   /* jump out of loop on error */

         for (ii = 0; ii < naxes[0]; ii++) {
           sum += pix[ii];                      /* accumlate sum */
           if (pix[ii] < minval) minval = pix[ii];  /* find min and  */
           if (pix[ii] > maxval) maxval = pix[ii];  /* max values    */
         }
      }
      
      free(pix);
      fits_close_file(fptr, &status);
    }

    if (status)	sciprint("operate error !"); /* print any error message */   
	else if (totpix > 0) meanval = sum / totpix;

	out_axes =(int*)malloc(2*sizeof(int));
	out_axes[0]=naxes[0]; out_axes[1]=naxes[1];
	CreateVarFromPtr(2, "i", &m2, &n2, &out_axes);   
	LhsVar(1) = 2;

	outptr=(double*)malloc(4*sizeof(double));
	outptr[0]=sum;outptr[1]=meanval;
	outptr[2]=minval;outptr[3]=maxval;
	CreateVarFromPtr(3, "d", &m3, &n3, &outptr);   
	LhsVar(2) = 3;
	return 0;
}


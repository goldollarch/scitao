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

#include "arroyo_wrap.h"
using namespace Arroyo;

RECONSTRUCTOR create_RECONSTRUCTOR ( APERTURE AP,
	DeformableMirror DM, LensletArray LenArray, ZERNIKE projected_modes, 
	double subaperture_illumination_threshold, double eigenvalue_threshold )
{
	RECONSTRUCTOR Reconstructor;

	Reconstructor.AP = AP; 
	Reconstructor.DM = DM;
	Reconstructor.LenArray = LenArray;
	Reconstructor.projected_modes = projected_modes;
	Reconstructor.eigenvalue_threshold=eigenvalue_threshold;
	Reconstructor.subaperture_illumination_threshold=subaperture_illumination_threshold;

	aperture *ap = c2cpp_APERTURE( AP );

	vector<long> actuator_dimensions(2);
	actuator_dimensions[0]=DM.actuator_x_dim; 
	actuator_dimensions[1]=DM.actuator_y_dim;
    long nactuators = ( DM.actuator_x_dim ) * ( DM.actuator_y_dim );

	square_lenslet_array temp_lnslt_arr	=c2cpp_LensletArray(LenArray);
	long nlenslets = ( LenArray.lenslet_x_axes ) * ( LenArray.lenslet_y_axes );

    vector<long> geometry_matrix_axes(2,2*nlenslets);
    geometry_matrix_axes[1] = nactuators;

	zernike znke=c2cpp_zernike( projected_modes );
	
	vector<long> reconstructor_axes(geometry_matrix_axes);
    reconstructor_axes[1] += znke.get_axes()[0];

	int reconstructor_matrix_nelem=reconstructor_axes[0]*reconstructor_axes[1];
	Reconstructor.reconstructor_matrix=create_PixelArray(reconstructor_axes[0],reconstructor_axes[1],0);
	
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	switch(DM.AP.aperture_type)  {
		case 1: 
			{
				in_diameter=DM.AP.ap_data2;
				out_diameter=DM.AP.ap_data1;
				ideal_deformable_mirror<annular_aperture>
					deformable_mirror (annular_aperture(in_diameter,out_diameter),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity);

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
		case 2:
			{
				x_size=DM.AP.ap_data1;
				y_size=DM.AP.ap_data2;
				ideal_deformable_mirror<rectangular_aperture> 
					deformable_mirror (rectangular_aperture(x_size,y_size),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity);

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
		case 3:
			{
				in_edge_length=DM.AP.ap_data1;
				ideal_deformable_mirror<hexagonal_aperture> 
					deformable_mirror (hexagonal_aperture(in_edge_length),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity);

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
		case 4:
			{
				in_diameter=DM.AP.ap_data2;
				out_diameter=DM.AP.ap_data1;
				nspiders=(int)(DM.AP.ap_data3); 
				spider_width=DM.AP.ap_data4;
				ideal_deformable_mirror<spidered_annular_aperture> 
					deformable_mirror (spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity);

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
/*
		case 5:
			{
  			    in_diameter=DM.AP.ap_data2;
			    out_diameter=DM.AP.ap_data1;
			    in_edge_length=DM.AP.ap_data3;	
			    in_gap_size=DM.AP.ap_data4;

				ideal_deformable_mirror<tiled_hexagonal_aperture> 
					deformable_mirror (tiled_hexagonal_aperture
					(in_diameter,out_diameter,in_edge_length,in_gap_size),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity);

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
*/
			//can not work! I don't know why!
		default:
			{
				diameter=DM.AP.ap_data1;
				ideal_deformable_mirror<circular_aperture> 	
					deformable_mirror ( circular_aperture( diameter ),
					actuator_dimensions, DM.actuator_pitch, DM.actuator_velocity );

				arroyo_least_squares_reconstructor<float> alsq_reconstructor
					(  *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					subaperture_illumination_threshold, eigenvalue_threshold, 0,0,0);

				for(int i=0;i<reconstructor_matrix_nelem;i++)
					Reconstructor.reconstructor_matrix.pixeldata[i]=alsq_reconstructor.data(i);

				break;
			}
	}

	return Reconstructor;
}

ZERNIKE arroyo_reconstruct_zernike_residuals( RECONSTRUCTOR Reconstructor,
					SHartmannCentroids SHcentroids)
{
	int xindex, yindex;
    double xresid, yresid;

	vector<long> actuator_axes(2);
	actuator_axes[0]=Reconstructor.DM.actuator_x_dim; 
	actuator_axes[1]=Reconstructor.DM.actuator_y_dim;

	vector<long> lenslet_axes(2);
	lenslet_axes[0] = Reconstructor.LenArray.lenslet_x_axes;
	lenslet_axes[1] = Reconstructor.LenArray.lenslet_y_axes;

	vector<long> centroid_axes( lenslet_axes );
	centroid_axes[1] *= 2;

    xindex = centroid_axes[0]*centroid_axes[1]*(actuator_axes[0]*actuator_axes[1]+1);
    yindex = xindex + centroid_axes[0]*centroid_axes[1];

	Shack_Hartmann_centroids shcentroids(c2cpp_SHartmannCentroids(SHcentroids));

    xresid = yresid = 0;
    for( int j=0; j< Reconstructor.reconstructor_matrix.x_axes; j++) {
		xresid -= Reconstructor.reconstructor_matrix.pixeldata[xindex+j]*shcentroids.data(j);
		yresid -= Reconstructor.reconstructor_matrix.pixeldata[yindex+j]*shcentroids.data(j);
    }

	zernike tip_tilt_residuals(1);	

    tip_tilt_residuals.set_cos_coeff(1,1,xresid);
    tip_tilt_residuals.set_sin_coeff(1,1,yresid);

	return(cpp2c_zernike(tip_tilt_residuals));
}

PixelArray arroyo_reconstruct_zonal_residuals( RECONSTRUCTOR Reconstructor,
					const SHartmannCentroids SHcentroids)
{
    int app_sign = 1;//    app_sign = app_sign_convention ? -1 : 1;
	long index = 0; double resid;

	vector<long> actuator_axes(2);
	actuator_axes[0]=Reconstructor.DM.actuator_x_dim; 
	actuator_axes[1]=Reconstructor.DM.actuator_y_dim;

	vector<long> lenslet_axes(2);
	lenslet_axes[0] = Reconstructor.LenArray.lenslet_x_axes;
	lenslet_axes[1] = Reconstructor.LenArray.lenslet_y_axes;

	long nsubaps = lenslet_axes[0]*lenslet_axes[1];
    int nactuators = actuator_axes[0]*actuator_axes[1];

	Shack_Hartmann_centroids shcentroids(c2cpp_SHartmannCentroids(SHcentroids));

	pixel_array<double> dm_residuals(actuator_axes);

    for(int i=0; i<actuator_axes[1]; i++)	{
		for(int j=0; j<actuator_axes[0]; j++)	{
			resid = 0;
			for(int k=0; k<nsubaps; k++)  {
                resid -= app_sign*Reconstructor.reconstructor_matrix.pixeldata[index+k]*shcentroids.data(k);
				resid -= Reconstructor.reconstructor_matrix.pixeldata[index+nsubaps+k]*shcentroids.data(k+nsubaps);
			}
			
			dm_residuals.set_data( j*actuator_axes[1]+i, resid );
			index += Reconstructor.reconstructor_matrix.x_axes;
		}
	}

	return(cpp2c_PixelArray(dm_residuals));
}

zernike_projected_zonal_reconstructor
*c2cpp_RECONSTRUCTOR( RECONSTRUCTOR Reconstructor)
{
	int nspiders;

	double diameter,in_diameter,out_diameter,
		x_size,y_size,in_edge_length,spider_width,in_gap_size;

	zernike_projected_zonal_reconstructor * zpz_reconstructor;

	aperture *ap = c2cpp_APERTURE( Reconstructor.AP );

	square_lenslet_array temp_lnslt_arr
		=c2cpp_LensletArray(Reconstructor.LenArray);

	zernike znke=c2cpp_zernike( Reconstructor.projected_modes );

	vector<long> actuator_dimensions(2);
	actuator_dimensions[0]=Reconstructor.DM.actuator_x_dim; 
	actuator_dimensions[1]=Reconstructor.DM.actuator_y_dim;
	
	switch(Reconstructor.DM.AP.aperture_type)  {
		case 1: 
			{
				in_diameter=Reconstructor.DM.AP.ap_data2;
				out_diameter=Reconstructor.DM.AP.ap_data1;
				ideal_deformable_mirror<annular_aperture>
					deformable_mirror (annular_aperture(in_diameter,out_diameter),
					actuator_dimensions, Reconstructor.DM.actuator_pitch,
					Reconstructor.DM.actuator_velocity);

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
		case 2:
			{
				x_size=Reconstructor.DM.AP.ap_data1;
				y_size=Reconstructor.DM.AP.ap_data2;
				ideal_deformable_mirror<rectangular_aperture> 
					deformable_mirror (rectangular_aperture(x_size,y_size),
					actuator_dimensions, Reconstructor.DM.actuator_pitch,
					Reconstructor.DM.actuator_velocity);

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
		case 3:
			{
				in_edge_length=Reconstructor.DM.AP.ap_data1;
				ideal_deformable_mirror<hexagonal_aperture> 
					deformable_mirror (hexagonal_aperture(in_edge_length),
					actuator_dimensions, Reconstructor.DM.actuator_pitch,
					Reconstructor.DM.actuator_velocity);

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
		case 4:
			{
				in_diameter=Reconstructor.DM.AP.ap_data2;
				out_diameter=Reconstructor.DM.AP.ap_data1;
				nspiders=(int)(Reconstructor.DM.AP.ap_data3); 
				spider_width=Reconstructor.DM.AP.ap_data4;
				ideal_deformable_mirror<spidered_annular_aperture> 
					deformable_mirror (spidered_annular_aperture
					(in_diameter,out_diameter,nspiders,spider_width),
					actuator_dimensions, Reconstructor.DM.actuator_pitch, 
					Reconstructor.DM.actuator_velocity);

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
			/*	case 5:
			{
  			    in_diameter=Reconstructor.DM.AP.ap_data2;
			    out_diameter=Reconstructor.DM.AP.ap_data1;
			    in_edge_length=Reconstructor.DM.AP.ap_data3;	
			    in_gap_size=Reconstructor.DM.AP.ap_data4;

				ideal_deformable_mirror<tiled_hexagonal_aperture> 
					deformable_mirror (tiled_hexagonal_aperture
					(in_diameter,out_diameter,in_edge_length,in_gap_size),
					actuator_dimensions, Reconstructor.DM.actuator_pitch, 
					Reconstructor.DM.actuator_velocity);

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
			*/	
			//can not work! I don't know why!
		default:
			{
				diameter=Reconstructor.DM.AP.ap_data1;
				ideal_deformable_mirror<circular_aperture> 	
					deformable_mirror ( circular_aperture( diameter ),
					actuator_dimensions, Reconstructor.DM.actuator_pitch,
					Reconstructor.DM.actuator_velocity );

				zpz_reconstructor=  new arroyo_least_squares_reconstructor<float>
					( *ap, deformable_mirror, temp_lnslt_arr, znke, false,	
					Reconstructor.subaperture_illumination_threshold, 
					Reconstructor.eigenvalue_threshold, 0,0,0);

				break;
			}
	}

	// only float supported ! double can't run! I don't know why!

	return zpz_reconstructor;
}

void write_RECONSTRUCTOR_file(RECONSTRUCTOR Reconstructor, char *fname)
{
	zernike_projected_zonal_reconstructor 
		*zpz_reconstructor = c2cpp_RECONSTRUCTOR( Reconstructor );
	
	stringstream filename_stream;
	filename_stream.str("");filename_stream << fname << ".fits";
    zpz_reconstructor->write(filename_stream.str().c_str());
}

arroyo_least_squares_reconstructor<float> 
arroyo_least_squares_reconstructor_file_constructor( char *fname )
{
	stringstream filename_stream;
	filename_stream.str(""); filename_stream << fname << ".fits";
	arroyo_least_squares_reconstructor<float>  
		reconstructor(filename_stream.str().c_str());
	//  it can not run under debug mode. I don't know why! 

	return reconstructor;
}

long total_Reconstructor_array_space( RECONSTRUCTOR Reconstructor )
{
	int AP_nelem,DM_nelem,LenArray_nelem,ZNK_nelem;
	long PixelArray_nelem,output;

    AP_nelem=17; 
	
	int actuator_x_dim = Reconstructor.DM.actuator_x_dim; 
	int actuator_y_dim = Reconstructor.DM.actuator_y_dim;
	int nactuator = actuator_x_dim*actuator_y_dim;
	DM_nelem=2*nactuator+22;

	LenArray_nelem = 6;
	ZNK_nelem = total_znk_space(Reconstructor.projected_modes.order)+1;

	PixelArray_nelem = ( Reconstructor.reconstructor_matrix.x_axes) 
	  * ( Reconstructor.reconstructor_matrix.y_axes ) + 2;

	output=AP_nelem+DM_nelem+LenArray_nelem+ZNK_nelem+PixelArray_nelem+2;
	return output;

}

long total_Reconstructor_array_nelem(double *inptr)
{
	long k;

	k=2;
	int AP_nelem=17; 
	k += AP_nelem;

	int DM_nelem = 2*((int)inptr[k])*((int)inptr[k+1])+22;
	k += DM_nelem;

	int LenArray_nelem = 6; 
	k += LenArray_nelem;

	int ZNK_nelem=total_znk_space((int)inptr[k])+1;
	k += ZNK_nelem;

	long PixelArray_nelem = ((int)inptr[k])*((int)inptr[k+1])+2; 
	k += PixelArray_nelem;

	return k;
}

double *RECONSTRUCTOR2array ( RECONSTRUCTOR Reconstructor )
{
	long i, k, nelem;
	double *tmp, *outptr;

	nelem=total_Reconstructor_array_space(Reconstructor);
	outptr=(double*)calloc(nelem,sizeof(double));

	outptr[0]=Reconstructor.eigenvalue_threshold;
	outptr[1]=Reconstructor.subaperture_illumination_threshold;

	k=2;
	tmp=APERTURE2array(Reconstructor.AP);
	for(i=0;i<17;i++) {
		outptr[k]=tmp[i];
		k++;
	}

	int actuator_x_dim = Reconstructor.DM.actuator_x_dim; 
	int actuator_y_dim = Reconstructor.DM.actuator_y_dim;
	int nactuator = actuator_x_dim*actuator_y_dim;
	int DM_nelem=2*nactuator+22;
	tmp=DeformableMirror2array(Reconstructor.DM);
	for(i=0;i<DM_nelem;i++) {
		outptr[k]=tmp[i];
		k++;
	}

	tmp=LensletArray2array(Reconstructor.LenArray);
	for(i=0;i<6;i++) {
		outptr[k]=tmp[i];
		k++;
	}

	int ZNK_nelem=total_znk_space(Reconstructor.projected_modes.order)+1;
	tmp=ZERNIKE2array(Reconstructor.projected_modes);
	for(i=0;i<ZNK_nelem;i++) {
		outptr[k]=tmp[i];
		k++;
	}

	long PixelArray_nelem = ( Reconstructor.reconstructor_matrix.x_axes) 
	  * ( Reconstructor.reconstructor_matrix.y_axes ) + 2;
	tmp=PixelArray2array(Reconstructor.reconstructor_matrix);
	for(i=0;i<PixelArray_nelem;i++) {
		outptr[k]=tmp[i];
		k++;
	}

	return outptr;
}

RECONSTRUCTOR array2RECONSTRUCTOR( double *inptr)
{
	long k;
	RECONSTRUCTOR Reconstructor;

	Reconstructor.eigenvalue_threshold=inptr[0];
	Reconstructor.subaperture_illumination_threshold=inptr[1];

	k=2;
	APERTURE AP = array2APERTURE( inptr+k );
	Reconstructor.AP=AP;
	int AP_nelem=17; 
	k += AP_nelem;

	DeformableMirror DM = array2DeformableMirror( inptr+k );
	Reconstructor.DM=DM;
	int DM_nelem = 2*DM.actuator_x_dim*DM.actuator_y_dim+22;
	k += DM_nelem;

	LensletArray LenArray=array2LensletArray(inptr+k);
	Reconstructor.LenArray=LenArray;
	k += 6;

	ZERNIKE projected_modes=array2ZERNIKE(inptr+k);
	Reconstructor.projected_modes=projected_modes;
	int ZNK_nelem=total_znk_space(projected_modes.order)+1;
	k += ZNK_nelem;

	PixelArray reconstructor_matrix=array2PixelArray(inptr+k);
	Reconstructor.reconstructor_matrix=reconstructor_matrix;

	return Reconstructor;
}


RECONSTRUCTOR RECONSTRUCTOR_file_constructor(char *fname)
{
	char *file_all_name;
	RECONSTRUCTOR Reconstructor;

    fitsfile *fptr;     
    int value;
    char oldvalue[FLEN_VALUE], comment[FLEN_COMMENT];
    int status = 0;  
    int iomode= READONLY, keytype;
	long rec_mat_x, rec_mat_y, totpix, fpixel[2];

    double *pix;

    fits_open_file(&fptr, fname, iomode, &status);

	fits_read_key(fptr,TINT,"NAXIS1", &value, comment, &status);
	rec_mat_x = value; 

	fits_read_key(fptr,TINT,"NAXIS2", &value, comment, &status);
	rec_mat_y = value;

	fits_read_key(fptr,TINT,"ACTUATR1", &value, comment, &status);
	Reconstructor.DM.actuator_x_dim = value; 

	fits_read_key(fptr,TINT,"ACTUATR2", &value, comment, &status);
	Reconstructor.DM.actuator_y_dim = value; 

	fits_read_key(fptr,TINT,"CNTROID1", &value, comment, &status);
	Reconstructor.LenArray.lenslet_x_axes = value; 

	fits_read_key(fptr,TINT,"CNTROID2", &value, comment, &status);
	Reconstructor.LenArray.lenslet_y_axes = value/2; 

	totpix=rec_mat_x*rec_mat_y;
	pix = (double *) malloc( totpix * sizeof(double)); /* memory */
	if (pix == NULL) exit(1);

	fpixel[0] = 1; fpixel[1] = 1;
	fits_read_pix(fptr, TDOUBLE, fpixel, totpix, 0, pix, 0, &status);

	fits_close_file(fptr, &status);

	PixelArray rec_mat =create_PixelArray(rec_mat_x, rec_mat_y,0);
	set_PixelArray_data(&rec_mat,pix);

	Reconstructor.reconstructor_matrix=rec_mat;

	return Reconstructor;
}

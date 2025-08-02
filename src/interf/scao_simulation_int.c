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

#include "time.h"
#include "int_optics.h"

int scao_simulation_int(char *fname) 
{
	clock_t ci;

	int m1, n1, l1, m2, n2, l2, i,j,k, num, num_detected_emitters, 
	  minlhs=1, maxlhs=1, minrhs=0, maxrhs=2;
	SIMULATION_PARAMETER scao_para;
	char *rec_fname=NULL;

	double default_para[57]={ 
		25, 5, 1, 0, 0.002,		1,
		100, 0, 5e-2, 2.67e-4, 0, 1e-7,
		0.6e-6, 1,		2, 3, 0, 0,
		10000, 5000, 5, 30,		0.02, 0.02,
		1, 1, 0,		1, 1,
		0, 0, 2, 0, 0, 0, 0,
		0.012246, 0.000252, 32, 32,
		0, 1,		16, 0.5, 1e-6,
		1, 0, 0.05, 0.5e-6, 0, 0, 0, 0, 0,
		2.2e-6, 8, 8
	};

	double guide_star_atmospheric_propagation_time=0;
	double science_target_atmospheric_propagation_time=0;
	double ttm_dm_correction_time=0;
	double far_field_propagation_time=0;
	double lenslet_array_transformation_time=0;
	double centroid_time=0;
	double reconstruction_time=0;
	double update_time=0;
	double t;

	CheckRhs(minrhs,maxrhs);
	CheckLhs(minlhs,maxlhs);

	scao_para=array2SimuPara(default_para);

	if (Rhs>0) {
		GetRhsVar(1, "d", &m1, &n1, &l1);
		num=54+((int)(*stk(l1+46)))+3*((int)(*stk(l1+13)));
		if ( (m1!=1)||(n1 != num ) )  {
			sciprint("Error: first arguments must be simulation parameter\r\n");
			return 0;
		}
		scao_para=array2SimuPara(stk(l1));
	}

	if (Rhs>1) {
		GetRhsVar(2, "c", &m2, &n2, &l2);
		rec_fname=cstk(l2);
	}
    
	scao_aperture(scao_para);

	if(scao_para.verbose) 
		sciprint("Constructing reconstructor...");
	ci = clock();
	scao_reconstructor(scao_para,rec_fname);
	if(scao_para.verbose) 
		sciprint("required %f seconds\n",(clock()-ci)/(double)CLOCKS_PER_SEC); 

	scao_ref_atm_model(scao_para);

	if(scao_para.verbose) 
		sciprint("Constructing guide star emitter and science emitters...\n");
	scao_emitter(scao_para);
	scao_wavefront_header(scao_para);

	if(scao_para.verbose) 
		sciprint("Constructing refractive atmospheric layers...");
	ci = clock();
	scao_ref_atm_layer(scao_para);
	if(scao_para.verbose) 
		sciprint("required %f seconds\n",(clock()-ci)/(double)CLOCKS_PER_SEC); 

	if(scao_para.verbose) 
		sciprint("Constructing tip tilt mirror and tip tilt proportional integral controller\n");
	scao_ttm_and_pi(scao_para);

	if(scao_para.verbose) 
		sciprint("Constructing deformable mirror and deformable proportional integral controller\n");
	scao_dm_and_pi(scao_para);

	if(scao_para.verbose) 
		sciprint("Constructing lenslet array\n");
	scao_lenslet_array(scao_para);


	if(scao_para.verbose)
		sciprint("Constructing unaberrated images\n");
	scao_unaberrated_image(scao_para);
	if(scao_para.verbose) 
		scao_write_wavefront(1,"unaberrated_focal_wf_", -1, 0);

	num_detected_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;

	for( i=0; i<scao_para.nsteps_in_simulation; i++)
	{
		sciprint("Timestep %d:\n",i);
		t=i*scao_para.timestep_seconds;

		if(i>=scao_para.nsteps_to_delay_propagating_science_wavefronts &&
			i%scao_para.interval_for_propagating_science_wavefronts==0)
		{
				for(j=0; j<scao_para.detected_wavelengths_number; j++)
				{
					for( k=0; k<num_detected_emitters; k++)
					{
						if(scao_para.verbose)
							sciprint("\tPropagating wavefront at wavelength %e through atmosphere...",
								scao_para.detected_wavelengths_meters[j]);
						ci = clock();
						scao_detected_aperture(scao_para,i,j,k);
						science_target_atmospheric_propagation_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
						if(scao_para.vverbose) 
							scao_write_wavefront(0,"uncorrected_pupil_wf_",t,k);
						if(scao_para.verbose) 
							sciprint("average time %f seconds\n",
								science_target_atmospheric_propagation_time/\
								(double)((i-scao_para.nsteps_to_delay_propagating_science_wavefronts+1)/\
								scao_para.interval_for_propagating_science_wavefronts));

						scao_uncorrected_detected_far_field(scao_para,i,j,k);
						scao_write_wavefront(1,"uncorrected_focal_wf_", t, k);

						if(scao_para.verbose)
							sciprint("\tCorrecting detected wavefront using tip tilt and deformable mirrors...\n");
						scao_correction_detected(scao_para,i,j,k);
						if(scao_para.verbose)
							scao_write_wavefront(0,"corrected_pupil_wf_", t, k);
						if(scao_para.vverbose) {
							if(k==0) {
								scao_write_wavefront(2,"tt_corrected_pupil_wf_", t, k);
								scao_write_wavefront(3,"dm_corrected_pupil_wf_", t, k);
							}
						}

						if(scao_para.verbose)
							sciprint("\tPropagating detected wavefront to far field...");
						scao_corrected_detected_far_field(scao_para,i,j,k);
						far_field_propagation_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
						if(scao_para.verbose)
							sciprint("average time %f seconds\n",far_field_propagation_time/(double)(i+1));
						scao_write_wavefront(0,"corrected_focal_wf_", t, k);
					}
				}
		}

		if(scao_para.verbose)
			sciprint("\tPropagating sensing wavefront through atmosphere...");
		ci = clock();
		scao_sensing_aperture(scao_para,i);
		guide_star_atmospheric_propagation_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) {
			sciprint("average time %f seconds\n",guide_star_atmospheric_propagation_time/(double)(i+1));
			scao_write_wavefront(0,"uncorrected_sensing_pupil_wf_", t, 0);
		}

		if(scao_para.verbose) 
			sciprint("\tCorrecting sensing wavefront using tip tilt and deformable mirrors...");
		ci = clock();    
		scao_correction_sensing(scao_para,i);
		ttm_dm_correction_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose)
			sciprint("average time %f seconds\n",ttm_dm_correction_time/(double)(i+1));
		if(scao_para.vverbose) {
			scao_write_wavefront(2,"tt_sensing_corrected_pupil_wf_", t, 0);
			scao_write_wavefront(3,"dm_sensing_corrected_pupil_wf_", t, 0);
		}

		if(scao_para.verbose) 
			sciprint("\tPropagating sensing wavefront through lenslet array...");
		ci = clock();
		scao_sensing_lenslet_array(scao_para,i);
		lenslet_array_transformation_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			sciprint("average time %f seconds\n",lenslet_array_transformation_time/(double)(i+1));
		if(scao_para.vverbose) 
			scao_write_wavefront(0,"lnslt_wf_", t, 0);

		if(scao_para.verbose) 
			sciprint("\tConstructing Shack Hartmann centroids...");
		ci = clock();
		scao_construct_centroids(scao_para,i);
		centroid_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			sciprint("average time %f seconds\n",centroid_time/(double)(i+1));

		if(scao_para.verbose) 
			sciprint("\tPerforming reconstruction...");
		ci = clock();
		scao_reconstruct_residuals(scao_para,i);
		reconstruction_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			sciprint("average time %f seconds\n",reconstruction_time/(double)(i+1));

		if(scao_para.verbose) 
			sciprint("\tUpdating tip tilt and deformable mirrors...");
		ci = clock();
		scao_update_corrector_mirrors(scao_para,i);
		update_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			sciprint("average time %f seconds\n",update_time/(double)(i+1));

		t= (guide_star_atmospheric_propagation_time + ttm_dm_correction_time +
				far_field_propagation_time+ lenslet_array_transformation_time + 
				centroid_time + reconstruction_time + update_time)/(double)(i+1);

		if(scao_para.verbose)
			sciprint("\tAverage time per timestep %f\n",t);
    
  }

  sciprint( "Average guide star atmospheric propagation time %f\n",
        guide_star_atmospheric_propagation_time/(double)scao_para.nsteps_in_simulation);

  sciprint( "Average tip tilt and deformable mirror correction time %f\n",
        ttm_dm_correction_time/(double)scao_para.nsteps_in_simulation );
  sciprint( "Average far field propagation time %f\n",
        far_field_propagation_time/(double)scao_para.nsteps_in_simulation );
  sciprint( "Average lenslet array transformation time %f\n",
        lenslet_array_transformation_time/(double)scao_para.nsteps_in_simulation );
  sciprint( "Average centroid computation time %f\n",
        centroid_time/(double)scao_para.nsteps_in_simulation );  
  sciprint( "Average reconstruction time %f\n",
        reconstruction_time/(double)scao_para.nsteps_in_simulation );
  sciprint( "Average mirror update time %f\n",
        update_time/(double)scao_para.nsteps_in_simulation );

  t = guide_star_atmospheric_propagation_time + ttm_dm_correction_time +
		far_field_propagation_time+ lenslet_array_transformation_time + 
		centroid_time + reconstruction_time + update_time;
  sciprint( "Average time per timestep %f\n",
			t/(double)scao_para.nsteps_in_simulation);
  
  return 0;

}

/*
SAO - A toolbox based on Scilab/Scicos environment for the simulation 
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

#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>

#include "wrapper/simulation.h"

clock_t ci;

void scao_initialize( SIMULATION_PARAMETER scao_para, char *rec_file_name )
{
	scao_aperture(scao_para);

	if(scao_para.verbose) 
		printf("Constructing reconstructor...");
	ci = clock();
	scao_reconstructor(scao_para,rec_file_name);
	if(scao_para.verbose) 
		printf("required %f seconds\n",(clock()-ci)/(double)CLOCKS_PER_SEC); 

	scao_ref_atm_model(scao_para);

	if(scao_para.verbose) 
		printf("Constructing guide star emitter and science emitters...\n");
	scao_emitter(scao_para);
	scao_wavefront_header(scao_para);

	if(scao_para.verbose) 
		printf("Constructing refractive atmospheric layers...");
	ci = clock();
	scao_ref_atm_layer(scao_para);
	if(scao_para.verbose) 
		printf("required %f seconds\n",(clock()-ci)/(double)CLOCKS_PER_SEC); 

	if(scao_para.verbose) 
		printf("Constructing tip tilt mirror and tip tilt proportional integral controller\n");
	scao_ttm_and_pi(scao_para);

	if(scao_para.verbose) 
		printf("Constructing deformable mirror and deformable proportional integral controller\n");
	scao_dm_and_pi(scao_para);

	if(scao_para.verbose) 
		printf("Constructing lenslet array\n");
	scao_lenslet_array(scao_para);

}

void scao_simulation( SIMULATION_PARAMETER scao_para)
{

	double guide_star_atmospheric_propagation_time=0;
	double science_target_atmospheric_propagation_time=0;
	double ttm_dm_correction_time=0;
	double far_field_propagation_time=0;
	double lenslet_array_transformation_time=0;
	double centroid_time=0;
	double reconstruction_time=0;
	double update_time=0;

	int i,j,k,num_detected_emitters;
	double t;

	if(scao_para.verbose)
		printf("Constructing unaberrated images\n");
	scao_unaberrated_image(scao_para);
	num_detected_emitters=scao_para.n_x_emitters*scao_para.n_y_emitters;

	for( i=0; i<scao_para.nsteps_in_simulation; i++){
		printf("Timestep %d:\n",i);
		if(i>=scao_para.nsteps_to_delay_propagating_science_wavefronts &&
			i%scao_para.interval_for_propagating_science_wavefronts==0)
		{
				for(j=0; j<scao_para.detected_wavelengths_number; j++)
				{
					for( k=0; k<num_detected_emitters; k++)
					{
						if(scao_para.verbose)
							printf("\tPropagating wavefront from emitter %d at wavelength %e through atmosphere...",
								k,scao_para.detected_wavelengths_meters[j]);
						ci = clock();
						scao_detected_aperture(scao_para,i,j,k);
						science_target_atmospheric_propagation_time +=
							(clock()-ci)/(double)CLOCKS_PER_SEC;
						t=science_target_atmospheric_propagation_time/
							(double)((i-scao_para.nsteps_to_delay_propagating_science_wavefronts+1)/
							scao_para.interval_for_propagating_science_wavefronts);
						if(scao_para.verbose) 
							printf("average time %e seconds\n",t);

						scao_uncorrected_detected_far_field(scao_para,i,j,k);
						t=i*scao_para.timestep_seconds;
						scao_write_wavefront(1,"uncorrected_focal_wf_", t, k);

						if(scao_para.verbose)
							printf("\tCorrecting detected wavefront using tip tilt and deformable mirrors...\n");
						scao_correction_detected(scao_para,i,j,k);

						if(scao_para.verbose)
							printf("\tPropagating detected wavefront to far field...");
						scao_corrected_detected_far_field(scao_para,i,j,k);
						far_field_propagation_time += 
							(clock()-ci)/(double)CLOCKS_PER_SEC;
						if(scao_para.verbose)
							printf("average time %e seconds\n",
								far_field_propagation_time/(double)(i+1));

						t=i*scao_para.timestep_seconds;
						scao_write_wavefront(0,"corrected_focal_wf_", t, k);
					}
				}
		}

		if(scao_para.verbose)
			printf("\tPropagating sensing wavefront through atmosphere...");
		ci = clock();
		scao_sensing_aperture(scao_para,i);
		guide_star_atmospheric_propagation_time += 
			(clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			printf("average time %e seconds\n",
				guide_star_atmospheric_propagation_time/(double)(i+1));

		if(scao_para.verbose) 
			printf("\tCorrecting sensing wavefront using tip tilt and deformable mirrors...");
		ci = clock();    
		scao_correction_sensing(scao_para,i);
		ttm_dm_correction_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose)
			printf("average time %e seconds\n",ttm_dm_correction_time/(double)(i+1));

		if(scao_para.verbose) 
			printf("\tPropagating sensing wavefront through lenslet array...");
		ci = clock();
		scao_sensing_lenslet_array(scao_para,i);
		lenslet_array_transformation_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			printf("average time %e seconds\n",
				lenslet_array_transformation_time/(double)(i+1));

		if(scao_para.verbose) 
			printf("\tConstructing Shack Hartmann centroids...");
		ci = clock();
		scao_construct_centroids(scao_para,i);
		centroid_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			printf("average time %e seconds\n",centroid_time/(double)(i+1));

		if(scao_para.verbose) 
			printf("\tPerforming reconstruction...");
		ci = clock();
		scao_reconstruct_residuals(scao_para,i);
		reconstruction_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			printf("average time %e seconds\n",reconstruction_time/(double)(i+1));

		if(scao_para.verbose) 
			printf("\tUpdating tip tilt and deformable mirrors...");
		ci = clock();
		scao_update_corrector_mirrors(scao_para,i);
		update_time += (clock()-ci)/(double)CLOCKS_PER_SEC;
		if(scao_para.verbose) 
			printf("average time %e seconds\n",update_time/(double)(i+1));

		t= (guide_star_atmospheric_propagation_time + ttm_dm_correction_time                   +
	       far_field_propagation_time+ lenslet_array_transformation_time + 
	       centroid_time + reconstruction_time + update_time) / (double)(i+1);

		if(scao_para.verbose)
			printf("\tAverage time per timestep %f\n",t);
    
  }

  printf( "Average guide star atmospheric propagation time %f\n",
        guide_star_atmospheric_propagation_time/
			(double)scao_para.nsteps_in_simulation);

  printf( "Average tip tilt and deformable mirror correction time %f\n",
        ttm_dm_correction_time/(double)scao_para.nsteps_in_simulation );
  printf( "Average far field propagation time %f\n",
        far_field_propagation_time/(double)scao_para.nsteps_in_simulation );
  printf( "Average lenslet array transformation time %f\n",
        lenslet_array_transformation_time/(double)scao_para.nsteps_in_simulation );
  printf( "Average centroid computation time %f\n",
        centroid_time/(double)scao_para.nsteps_in_simulation );  
  printf( "Average reconstruction time %f\n",
        reconstruction_time/(double)scao_para.nsteps_in_simulation );
  printf( "Average mirror update time %f\n",
        update_time/(double)scao_para.nsteps_in_simulation );

  t = t / (double) scao_para.nsteps_in_simulation;
  printf( "Average time per timestep %f\n",t);

}

void main()
{
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

	SIMULATION_PARAMETER para;
	para=array2SimuPara(default_para);

	scao_initialize(para,NULL);
	scao_simulation(para);

}

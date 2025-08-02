##
## -------------------------------------------------------------------------
## scitao - Scilab/Scicos Adaptive Optics toolbox
## Copyright (C) 2006  IAPCM, Beijing, China.
## Some inspiration came from SIP's contribution 
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  
## See the GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software Foundation,
## Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
## -------------------------------------------------------------------------
##

libdir=@TOOLBOXDIR@/lib
lib_LTLIBRARIES = libarroyo.la
EXTRA_DIST = builder.sce
AM_CPPFLAGS = @LIBRARY_EXTRA_CPPFLAGS@

libarroyo_la_SOURCES = AO_algo.cpp AO_algo.h AO_cpp.h AO_sim_base.cpp AO_sim_base.h ao_simulation.h aperture.cpp aperture.h aplusplus_reconstructor.cpp aplusplus_reconstructor.h arcadia_mcao_reconstructor.cpp arcadia_mcao_reconstructor.h arroyo.h arroyo_least_squares_reconstructor.cpp arroyo_least_squares_reconstructor.h atmospheric_model.h Cn2_distribution.cpp Cn2_distribution.h colormap.cpp colormap.h computational_geometry.cpp computational_geometry.h conic_mirror.cpp conic_mirror.h covariance.cpp covariance.h deformable_mirror.cpp deformable_mirror.h Dereferenced.h detector.h diffractive_wavefront.cpp diffractive_wavefront.h Ellerbroek_Cerro_Pachon_model.cpp Ellerbroek_Cerro_Pachon_model.h Ellerbroek_Mauna_Kea_model.cpp Ellerbroek_Mauna_Kea_model.h emitter.cpp emitter.h fft_manager.cpp fft_manager.h fits_factory.h fits_header_data.cpp fits_header_data.h Gemini_GLAO_study_model.cpp Gemini_GLAO_study_model.h getopt_ref_atm_model.h Hufnagel_Valley_model.cpp Hufnagel_Valley_model.h iofits.cpp iofits.h lenslet_array.cpp lenslet_array.h linear_algebra.cpp linear_algebra.h modal_expansion.cpp modal_expansion.h observation.cpp observation.h optical_system.h optic.cpp optic.h output_plan.h PALAO_reconstructor.cpp PALAO_reconstructor.h Palomar_model.cpp Palomar_model.h pixel_amp_array.h pixel_array.cpp pixel_array.h pixel_phase_array.h power_spectrum.cpp power_spectrum.h propagation_plan.cpp propagation_plan.h proportional_integral_controller.h refractive_atmosphere.cpp refractive_atmosphere.h refractive_atmospheric_layer.cpp refractive_atmospheric_layer.h region_base.cpp region_base.h Shack_Hartmann_centroids.cpp Shack_Hartmann_centroids.h shack_hartmann.h sim_utils.cpp sim_utils.h SLCSAT_model.cpp SLCSAT_model.h special_functions.cpp special_functions.h structure_function.cpp structure_function.h subharmonic_method.cpp subharmonic_method.h three_frame.cpp three_frame.h three_point.cpp three_point.h three_transformation.cpp three_transformation.h three_vector.cpp three_vector.h tip_tilt_mirror.cpp tip_tilt_mirror.h TMT_SRD_v13_Cn2_model.cpp TMT_SRD_v13_Cn2_model.h user_defaults.h wavefront.cpp wavefront.h wavefront_header.cpp wavefront_header.h wavefront_phase_estimate.cpp wavefront_phase_estimate.h wavefront_sensor.h wind_model.cpp wind_model.h zernike.cpp zernike.h zernike_projected_zonal_reconstructor.cpp zernike_projected_zonal_reconstructor.h

libarroyo_la_LDFLAGS = -no-undefined -version-info 0:0:0 -release @PACKAGE_VERSION@

libarroyo_la_LIBADD = @LIBRARY_EXTRA_LIBS@
SUBDIRS = win-util


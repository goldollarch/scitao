
include_directories(
  ${CFITSIO_INCLUDE_DIR}
  ${FFTW_INCLUDE_DIR}
  ${LAPACK_INCLUDE_DIR})

set(ARROYO_HDRS 
${CMAKE_CURRENT_SOURCE_DIR}/AO_algo.h
${CMAKE_CURRENT_SOURCE_DIR}/AO_cpp.h
${CMAKE_CURRENT_SOURCE_DIR}/AO_sim_base.h
${CMAKE_CURRENT_SOURCE_DIR}/ao_simulation.h
${CMAKE_CURRENT_SOURCE_DIR}/aperture.h
${CMAKE_CURRENT_SOURCE_DIR}/aplusplus_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/arcadia_mcao_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/arroyo.h
${CMAKE_CURRENT_SOURCE_DIR}/arroyo_least_squares_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/atmospheric_model.h
${CMAKE_CURRENT_SOURCE_DIR}/Cn2_distribution.h
${CMAKE_CURRENT_SOURCE_DIR}/colormap.h
${CMAKE_CURRENT_SOURCE_DIR}/computational_geometry.h
${CMAKE_CURRENT_SOURCE_DIR}/conic_mirror.h
${CMAKE_CURRENT_SOURCE_DIR}/conic_section.h
${CMAKE_CURRENT_SOURCE_DIR}/covariance.h
${CMAKE_CURRENT_SOURCE_DIR}/deformable_mirror.h
${CMAKE_CURRENT_SOURCE_DIR}/Dereferenced.h
${CMAKE_CURRENT_SOURCE_DIR}/detector.h
${CMAKE_CURRENT_SOURCE_DIR}/diffractive_wavefront.h
${CMAKE_CURRENT_SOURCE_DIR}/Ellerbroek_Cerro_Pachon_model.h
${CMAKE_CURRENT_SOURCE_DIR}/Ellerbroek_Mauna_Kea_model.h
${CMAKE_CURRENT_SOURCE_DIR}/emitter.h
${CMAKE_CURRENT_SOURCE_DIR}/fft_manager.h
${CMAKE_CURRENT_SOURCE_DIR}/fits_factory.h
${CMAKE_CURRENT_SOURCE_DIR}/fits_header_data.h
${CMAKE_CURRENT_SOURCE_DIR}/Gemini_GLAO_study_model.h
${CMAKE_CURRENT_SOURCE_DIR}/getopt_ref_atm_model.h
${CMAKE_CURRENT_SOURCE_DIR}/Hufnagel_Valley_model.h
${CMAKE_CURRENT_SOURCE_DIR}/iofits.h
${CMAKE_CURRENT_SOURCE_DIR}/lenslet_array.h
${CMAKE_CURRENT_SOURCE_DIR}/linear_algebra.h
${CMAKE_CURRENT_SOURCE_DIR}/modal_expansion.h
${CMAKE_CURRENT_SOURCE_DIR}/NGAO_system_design_model.h
${CMAKE_CURRENT_SOURCE_DIR}/observation.h
${CMAKE_CURRENT_SOURCE_DIR}/optical_system.h
${CMAKE_CURRENT_SOURCE_DIR}/optic.h
${CMAKE_CURRENT_SOURCE_DIR}/output_plan.h
${CMAKE_CURRENT_SOURCE_DIR}/PALAO_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/Palomar_model.h
${CMAKE_CURRENT_SOURCE_DIR}/pixel_amp_array.h
${CMAKE_CURRENT_SOURCE_DIR}/pixel_array.h
${CMAKE_CURRENT_SOURCE_DIR}/pixel_phase_array.h
${CMAKE_CURRENT_SOURCE_DIR}/power_spectrum.h
${CMAKE_CURRENT_SOURCE_DIR}/profile_timeseries.h
${CMAKE_CURRENT_SOURCE_DIR}/propagation_plan.h
${CMAKE_CURRENT_SOURCE_DIR}/proportional_integral_controller.h
${CMAKE_CURRENT_SOURCE_DIR}/refractive_atmosphere.h
${CMAKE_CURRENT_SOURCE_DIR}/refractive_atmospheric_layer.h
${CMAKE_CURRENT_SOURCE_DIR}/region_base.h
${CMAKE_CURRENT_SOURCE_DIR}/Shack_Hartmann_centroids.h
${CMAKE_CURRENT_SOURCE_DIR}/shack_hartmann.h
${CMAKE_CURRENT_SOURCE_DIR}/sim_utils.h
${CMAKE_CURRENT_SOURCE_DIR}/SLCSAT_model.h
${CMAKE_CURRENT_SOURCE_DIR}/special_functions.h
${CMAKE_CURRENT_SOURCE_DIR}/structure_function.h
${CMAKE_CURRENT_SOURCE_DIR}/subharmonic_method.h
${CMAKE_CURRENT_SOURCE_DIR}/svd.h
${CMAKE_CURRENT_SOURCE_DIR}/three_frame.h
${CMAKE_CURRENT_SOURCE_DIR}/three_point.h
${CMAKE_CURRENT_SOURCE_DIR}/three_transformation.h
${CMAKE_CURRENT_SOURCE_DIR}/three_vector.h
${CMAKE_CURRENT_SOURCE_DIR}/time_val.h
${CMAKE_CURRENT_SOURCE_DIR}/tip_tilt_mirror.h
${CMAKE_CURRENT_SOURCE_DIR}/TMT_SRD_v13_Cn2_model.h
${CMAKE_CURRENT_SOURCE_DIR}/tomographic_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/user_defaults.h
${CMAKE_CURRENT_SOURCE_DIR}/wavefront.h
${CMAKE_CURRENT_SOURCE_DIR}/wavefront_header.h
${CMAKE_CURRENT_SOURCE_DIR}/wavefront_phase_estimate.h
${CMAKE_CURRENT_SOURCE_DIR}/wavefront_sensor.h
${CMAKE_CURRENT_SOURCE_DIR}/wind_model.h
${CMAKE_CURRENT_SOURCE_DIR}/zernike.h
${CMAKE_CURRENT_SOURCE_DIR}/zernike_projected_zonal_reconstructor.h
${CMAKE_CURRENT_SOURCE_DIR}/geometric_wavefront.h
)

set(ARROYO_SRCS  
${CMAKE_CURRENT_SOURCE_DIR}/AO_algo.cpp
${CMAKE_CURRENT_SOURCE_DIR}/AO_sim_base.cpp
${CMAKE_CURRENT_SOURCE_DIR}/aperture.cpp
${CMAKE_CURRENT_SOURCE_DIR}/aplusplus_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/arcadia_mcao_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/arroyo_least_squares_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Cn2_distribution.cpp
${CMAKE_CURRENT_SOURCE_DIR}/colormap.cpp
${CMAKE_CURRENT_SOURCE_DIR}/computational_geometry.cpp
${CMAKE_CURRENT_SOURCE_DIR}/conic_mirror.cpp
${CMAKE_CURRENT_SOURCE_DIR}/conic_section.cpp
${CMAKE_CURRENT_SOURCE_DIR}/covariance.cpp
${CMAKE_CURRENT_SOURCE_DIR}/deformable_mirror.cpp
${CMAKE_CURRENT_SOURCE_DIR}/diffractive_wavefront.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Ellerbroek_Cerro_Pachon_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Ellerbroek_Mauna_Kea_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/emitter.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fft_manager.cpp
${CMAKE_CURRENT_SOURCE_DIR}/fits_header_data.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Gemini_GLAO_study_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/getopt_ref_atm_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Hufnagel_Valley_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/iofits.cpp
${CMAKE_CURRENT_SOURCE_DIR}/lenslet_array.cpp
${CMAKE_CURRENT_SOURCE_DIR}/linear_algebra.cpp
${CMAKE_CURRENT_SOURCE_DIR}/modal_expansion.cpp
${CMAKE_CURRENT_SOURCE_DIR}/NGAO_system_design_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/observation.cpp
${CMAKE_CURRENT_SOURCE_DIR}/optic.cpp
${CMAKE_CURRENT_SOURCE_DIR}/PALAO_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Palomar_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/pixel_array.cpp
${CMAKE_CURRENT_SOURCE_DIR}/power_spectrum.cpp
${CMAKE_CURRENT_SOURCE_DIR}/profile_timeseries.cpp
${CMAKE_CURRENT_SOURCE_DIR}/propagation_plan.cpp
${CMAKE_CURRENT_SOURCE_DIR}/refractive_atmosphere.cpp
${CMAKE_CURRENT_SOURCE_DIR}/refractive_atmospheric_layer.cpp
${CMAKE_CURRENT_SOURCE_DIR}/region_base.cpp
${CMAKE_CURRENT_SOURCE_DIR}/Shack_Hartmann_centroids.cpp
${CMAKE_CURRENT_SOURCE_DIR}/sim_utils.cpp
${CMAKE_CURRENT_SOURCE_DIR}/SLCSAT_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/special_functions.cpp
${CMAKE_CURRENT_SOURCE_DIR}/structure_function.cpp
${CMAKE_CURRENT_SOURCE_DIR}/subharmonic_method.cpp
${CMAKE_CURRENT_SOURCE_DIR}/three_frame.cpp
${CMAKE_CURRENT_SOURCE_DIR}/three_point.cpp
${CMAKE_CURRENT_SOURCE_DIR}/three_transformation.cpp
${CMAKE_CURRENT_SOURCE_DIR}/three_vector.cpp
${CMAKE_CURRENT_SOURCE_DIR}/time_val.cpp
${CMAKE_CURRENT_SOURCE_DIR}/tip_tilt_mirror.cpp
${CMAKE_CURRENT_SOURCE_DIR}/TMT_SRD_v13_Cn2_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/tomographic_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/wavefront.cpp
${CMAKE_CURRENT_SOURCE_DIR}/wavefront_header.cpp
${CMAKE_CURRENT_SOURCE_DIR}/wavefront_phase_estimate.cpp
${CMAKE_CURRENT_SOURCE_DIR}/wind_model.cpp
${CMAKE_CURRENT_SOURCE_DIR}/zernike.cpp
${CMAKE_CURRENT_SOURCE_DIR}/zernike_projected_zonal_reconstructor.cpp
${CMAKE_CURRENT_SOURCE_DIR}/geometric_wavefront.cpp
)

install(FILES ${ARROYO_HDRS}
	      DESTINATION include/arroyo)

add_library(arroyo ${ARROYO_SRCS} ${ARROYO_HDRS})

set(libs tosc_lite
  ${CFITSIO_LIBRARY}
  ${FFTW_LIBRARY}
  ${LAPACK_LIBRARIES}
  )

Install(TARGETS arroyo
                RUNTIME DESTINATION bin/${CMAKE_BUILD_TYPE_INIT}  COMPONENT Runtime
                LIBRARY DESTINATION lib/${CMAKE_BUILD_TYPE_INIT}  COMPONENT Runtime
                ARCHIVE DESTINATION lib/${CMAKE_BUILD_TYPE_INIT}  COMPONENT Development)


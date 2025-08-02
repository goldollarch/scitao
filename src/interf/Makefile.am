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
lib_LTLIBRARIES = libinterf.la

AM_CPPFLAGS = @LIBRARY_EXTRA_CPPFLAGS@ -I../lightPipes -I../wrapper

libinterf_la_SOURCES = aperture_int.c atmospheric_model_int.c controller_int.c deformable_mirror_int.c emitter_int.c geometry_int.c int_arroyo.h int_lightPipes.h int_optics.c int_optics.h int_others.h lightPipes_int.c others_int.c propagation_int.c reconstructor_int.c refractive_layer_int.c scao_simulation_int.c shack_hartmann_int.c some_util_int.c tip_tilt_mirror_int.c wavefront_int.c wind_model_int.c zernike_int.c

libinterf_la_LDFLAGS = -no-undefined -version-info 0:0:0 -release @PACKAGE_VERSION@

libinterf_la_LIBADD = @LIBRARY_EXTRA_LIBS@

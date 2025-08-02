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
lib_LTLIBRARIES = libarroyo_wrap.la
EXTRA_DIST = builder.sce
AM_CPPFLAGS = @LIBRARY_EXTRA_CPPFLAGS@ -I../lightPipes -I../arroyo

libarroyo_wrap_la_SOURCES = aperture.cpp aperture.h arroyo_wrap.h atmospheric_model.cpp atmospheric_model.h controller.cpp controller.h deformable_mirror.cpp deformable_mirror.h emitter.cpp emitter.h geometry.cpp geometry.h propagation.cpp propagation.h reconstructor.cpp reconstructor.h refractive_layer.cpp refractive_layer.h scao_simulation.cpp shack_hartmann.cpp shack_hartmann.h simulation.cpp simulation.h some_utils.cpp some_utils.h tip_tilt_mirror.cpp tip_tilt_mirror.h wavefront.cpp wavefront.h wind_model.cpp wind_model.h zernike.cpp zernike.h

libarroyo_wrap_la_LDFLAGS = -no-undefined -version-info 0:0:0 -release @PACKAGE_VERSION@

libarroyo_wrap_la_LIBADD = @LIBRARY_EXTRA_LIBS@

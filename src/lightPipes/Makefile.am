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
lib_LTLIBRARIES = liblightpipes.la
EXTRA_DIST = builder.sce
AM_CPPFLAGS = @LIBRARY_EXTRA_CPPFLAGS@

liblightpipes_la_SOURCES =	 absorber.c begin.c b_mix.c b_split.c circ_ap.c circ_screen.c convert.c cros_out.c c_scilab.c fftn.c fftn.h fft_prop.c file_int.c file_pgm.c file_pha.c file_ps.c fil_ter.c forvard.c forward.c fresnel.c frsn_prop.c gauss.c gauss_screen.c interp1.c interpol.c int_pha.c l_amplif.c lens.c lens_forvard.c lens_fresn.c lightPipes.h normal.c pipes.c pipes.h pip_fft.c random.c rect_ap.c rect_screen.c steps.c strehl.c tilt.c tor_lens.c unfolds.c zernike.c

liblightpipes_la_LDFLAGS = -no-undefined -version-info 0:0:0 -release @PACKAGE_VERSION@

liblightpipes_la_LIBADD = @LIBRARY_EXTRA_LIBS@

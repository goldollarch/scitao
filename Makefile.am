##
## -------------------------------------------------------------------------
## scitao - Scilab/Scicos Adaptive Optics toolbox
## Copyright (C) 2006  IAPCM, Beijing, China.
## Some inspiration came from sip's contribution 
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

AUTOMAKE_OPTIONS = foreign
ACLOCAL_AMFLAGS = -I config

EXTRA_DIST = win-builder.sce builder.sce loader.sce \
	AUTHORS.txt COPYING.txt INSTALL.txt README.txt 

SUBDIRS = src macros man demos examples

toolboxdir=@TOOLBOXDIR@
toolbox_DATA = builder.sce loader.sce

distclean-local:
	rm -f path.incl

autoload:
	@echo "Appending code to $(HOME)/.scilab"
	@echo "//-- BEGIN CODE TO LOAD TOOLBOX --" >> $(HOME)/.scilab
	@echo "exec('@TOOLBOXDIR@/loader.sce');"   >> $(HOME)/.scilab
	@echo "//--------------------------------" >> $(HOME)/.scilab

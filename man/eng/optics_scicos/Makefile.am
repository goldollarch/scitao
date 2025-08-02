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

images  = absorb_em.gif absorber.gif animate.gif anim_em.gif ap_em.gif ap.gif atm_em.gif atm.gif clone.gif convert.gif diagnose_em.gif diagnose.gif dis_em.gif display.gif dm.gif emitter.gif emt_em.gif fft_em.gif fft.gif fits_em.gif fits.gif gauss_em.gif gauss.gif geo_em.gif geom.gif hartmann.gif interpolate.gif len_em.gif len.gif lens_em.gif mix_em.gif mixer.gif normal.gif norm_em.gif pic.gif prop.gif recon.gif resize_em.gif scao_em.gif split_em.gif split.gif start.gif ttm.gif znike.gif znk_em.gif

manpages_xml  = absorber.xml animate.xml aperture.xml atmosphere.xml clone.xml controller.xml convert.xml df_mirror.xml diagnose.xml display.xml emitter.xml fft.xml file_fits.xml gauss.xml geometry.xml hartmann.xml interpolate.xml lens.xml mixer.xml normal.xml propagate.xml reconstructor.xml splitter.xml start.xml tt_mirror.xml zernike.xml

manpages_html = $(manpages_xml:.xml=.htm)

helpdir=@TOOLBOXDIR@/man/eng/optics_scicos

help_DATA = $(manpages_html) $(images) whatis.htm
EXTRA_DIST = $(manpages_xml) $(images) scicos_rev.xsl

.xml.htm: $(manpages_xml)
	@echo "Creating html help files"
	@$(SCIDIR)/bin/scilab -nwni -nb -e "exec('../../sao_xmltohtml.sci');%helps=[PWD,'Optics_Scicos';%helps];sao_xmltohtml(PWD,'Optics_scicos','scicos_rev.xsl');quit"

clean-local:
	rm -f *.htm


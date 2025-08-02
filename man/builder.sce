mode(-1)
//
// -------------------------------------------------------------------------
// SAO - Scilab/Scicos Adaptive Optics tooolbox
//
// Copyright (C) 2006  IAPCM , Beijing, China.  Written by
// Chen jingyuan.  For comments or questions about this software,
// please contact the author at jingyuan_chen@yahoo.com.cn.
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software Foundation, 
// Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
// -------------------------------------------------------------------------
//

path=get_absolute_file_path('builder.sce');

myhelps=[path+'eng/optics', "Optics functions";
        path+'eng/optics_scicos', "Optics Scicos-blocks"]

exec('sao_xmltohtml.sci');

// update %helps for cross reference
%helps_save=%helps
%helps=[%helps;myhelps] 

if MSDOS then
	sao_xmltohtml(myhelps(1,1),myhelps(1,2),"scilab_rev.xsl");
	sao_xmltohtml(myhelps(2,1),myhelps(2,2),"scicos_rev.xsl");
end

sao_xmltohtml(myhelps(:,1),myhelps(:,2),"index");

//restore the previous help table
%helps=%helps_save
//erase temporary variables
clear %helps_save xmltohtml myhelps path get_absolute_file_path

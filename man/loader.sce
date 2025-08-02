mode(-1) //force silent execution
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

path=get_absolute_file_path('loader.sce');
//get the absolute path of this file

add_help_chapter("Optics functions",path+'eng/optics');
//add first help chapter
add_help_chapter("Optics Scicos-blocks",path+'eng/optics_scicos');
//add second help chapter

//clear the variable stack
clear path add_help_chapter get_absolute_file_path



	
	


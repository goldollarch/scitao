//
// -------------------------------------------------------------------------
// scitao - Scilab/Scicos Adaptive Optics tooolbox
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

function [] = field_plot(field)

  number=field(1); 
  x=1:number;
  
  number2=field(11); 
  y=1:number2;
  
  intensity=field_int(field);
  phase=field_pha(field);
 
  xbasc();
  drawlater();

  ncolor=60;

  subplot(2,2,1);
  h=hotcolormap(ncolor);
  xset('colormap',h);
  colorbar(0,1);
  grayplot(x,y,intensity,strf="030",rect=[1,1,number,number2]);

  subplot(2,2,2);
  h=hotcolormap(ncolor);
  xset('colormap',h);
  surf(intensity);
  xtitle( 'the intensty of the field', 'X', 'Y', 'intensity') ;

  subplot(2,2,3);
  h=hotcolormap(ncolor);
  xset('colormap',h);
  colorbar(0,1);
  grayplot(x,y,phase,strf="030",rect=[1,1,number,number2]);
  
  subplot(2,2,4);
  h=hotcolormap(ncolor);
  xset('colormap',h);
  surf(phase);
  xtitle( 'the phase of the field', 'X', 'Y', 'phase') ;

  drawnow();
  
endfunction 


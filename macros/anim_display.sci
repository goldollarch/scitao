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

function [blocks] = anim_display(blocks,flag)

  wid=10+curblock();
  
  if flag<>4 then
    graph=scf(wid);
  end
  
  ipar=blocks.ipar
  nipar=blocks.nipar
  
  iwp=nipar-3
  iwd=nipar-1
  
  num=ipar(1)
  wpx=ipar(2); 
  wpy=ipar(3);
  
  width=ipar(4); 
  height=ipar(5);
  
  x = 1:num;
  y = 1:num;
  z = zeros(num,num);

  if flag==4 then 
  
    xset("window",wid)
    if ipar(iwp) >=0 then 
       xset("wpos",wpx,wpy);
    end
    if ipar(iwd) >=0 then 
       xset("wpdim",width,height);
 //      xset("window",wid);
    end
    
    set("figure_style","new")

    graph=scf(wid); 
    clf(graph)
    
    graph.pixmap='on'
    Axe=graph.children
    
    Axe.isoview='on'
    
    tl=Axe.title;  
    tl.foreground=9; 
    tl.font_size=5; 
    tl.font_style=5;
    tl.text="field intensity";
//    tl.fill_mode="on"

    xset("colormap",hotcolormap(64))
//    xset("colormap",jetcolormap(64))
    
    colorbar(0,1)
    Sgrayplot(x,y,z,strf="042",zminmax=[0,1])

  elseif flag==2 then
  
    Axe=graph.children
    Axe.box='off'
    
    c=gce()
    e=c.children
    
    re=[];im=[]
    nelem=num*num
    for i=1:nelem
       re=[re,blocks.inptr(1)(i+30)];
       im=[im,blocks.inptr(1)(i+30+nelem)];
    end
    in=re^2+im^2
    e.data(:,3)=in'

    show_pixmap()

  end
  
endfunction 


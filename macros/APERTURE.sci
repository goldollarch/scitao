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

function [x,y,typ]=APERTURE(job,arg1,arg2)
// Copyright INRIA
x=[];y=[];typ=[];
select job
case 'plot' then
  standard_draw(arg1)
case 'getinputs' then
  [x,y,typ]=standard_inputs(arg1)
case 'getoutputs' then
  [x,y,typ]=standard_outputs(arg1)
case 'getorigin' then
  [x,y]=standard_origin(arg1)
case 'set' then
  x=arg1;
  graphics=arg1.graphics;exprs=graphics.exprs
  model=arg1.model;
  while %t do
    [ok,nport,ap_typ,data1,data2,data3,data4,exprs]=getvalue(..
       ['Set the aperture parameter'],..
	      ['the number of ports ( 1/2 )';
	      'aperture type ( 0/1/2/3/4/5 )';
	      'diameter,x_size,or in_edge_length';
             'in_diameter, or y_size';
             'nspiders, or in_edge_length';
             'spider_width, or in_gap_size'],..
	    list('vec',1,'vec',1,'vec',1,'vec',1,'vec',1,'vec',1),exprs)
    if ~ok then break,end
    if data1<0 then
      message('diameter,x_size,or in_edge_length must be positive')
      ok=%f
    end
    if data2<0 then
      message('in_diameter, or y_size must be positive')
      ok=%f
    end
    if data3<0 then
      message('nspiders, or in_edge_length must be positive')
      ok=%f
    end
    if data4<0 then
      message('spider_width, or in_gap_size must be positive')
      ok=%f
    end
    
    if nport==2 then
      [model,graphics,ok]=check_io(model,graphics,[-1;-1],[-1;-1],[],[])
    else
      [model,graphics,ok]=check_io(model,graphics,-1,-1,[],[])
    end
    
    if ok then
      ipar=ap_typ;
      rpar=[data1;data2;data3;data4];
      model.ipar=ipar;
      model.rpar=rpar;
      graphics.exprs=exprs;
      x.graphics=graphics;x.model=model
      break
    end
  end
case 'define' then
  model=scicos_model()
  model.sim=list('aperture',4)
  model.in=-1
  model.out=-1
  model.rpar=[5;2;4;0.2]
  model.ipar=0
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(nport);string(ap_typ);
       string(data1);string(data2);
       string(data3);string(data4)]
  gr_i=['txt=[''AP''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([3 2],model,exprs,gr_i)
end
endfunction

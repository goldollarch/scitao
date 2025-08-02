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

function [x,y,typ]=CLONE(job,arg1,arg2)
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
    [ok,in,exprs]=getvalue(['Set module parameter'],..
	    ['output ports number'],..
	    list('vec',1),exprs)
    if ~ok then break,end
   [model,graphics,ok]=check_io(model,graphics,-1,-ones(1,in),[],[])
    if ok then
      ipar=in;
      model.ipar=ipar;
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
case 'define' then
  model=scicos_model()
  model.sim=list('mixer',4)
  model.in=-1
  model.out=[-1;-1]
  model.ipar=2
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(in)]
  gr_i=''
  x=standard_define([.5 2],model,exprs,gr_i)
end
endfunction

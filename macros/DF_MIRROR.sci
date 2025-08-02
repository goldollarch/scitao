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

function [x,y,typ]=DF_MIRROR(job,arg1,arg2)
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
    [ok,ap,dim,pitch,vel,exprs]=getvalue(..
       ['Set the deformable mirror parameter'],..
	      ['mirror aperture parameter';
	      'the dimension of actuators';
	      'the actuator pitch (meter)';
	      'actuator velocity (m/s)'],..
	    list('vec',-1,'vec',-1,'vec',1,'vec',1),exprs)
	    
	if ~ok then break,end
    [model,graphics,ok]=check_io(model,graphics,[-1;-1;-1],[-1;-1],[],[])
    
    if ok then
      rpar=[ap(2);ap(3);ap(4);ap(5);pitch;vel];
      ipar=[ap(1);dim(:)]
      model.rpar=rpar;
      model.ipar=ipar;
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
case 'define' then
  model=scicos_model()
  model.sim=list('df_mirror',4)
  model.in=[-1;-1;-1]
  model.out=[-1;-1]
  model.rpar=[5;2;4;0.2;1;100]
  model.ipar=[0;8;8]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(ap);string(dim);string(pitch);string(vel)]
  gr_i=['txt=[''DM''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([3 2],model,exprs,gr_i)
end
endfunction

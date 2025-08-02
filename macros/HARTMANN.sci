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

function [x,y,typ]=HARTMANN(job,arg1,arg2)
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
    [ok,axes,pitch,focus,finpix,fname,exprs]=getvalue(..
       ['Set Shack_Hartmann sensor parameter'],..
	        ['the number of lenslets ';
	        'size of the square lenslet (meter)';
	        'focal lenth of the lenslets (meter)';
	        'final wavefront pixels per lenslet/transform ';
	        'transform wavefront file name'],..
	    list('vec',-1,'vec',1,'vec',1,'vec',-1,'str',1),exprs)
    if ~ok then break,end
    [model,graphics,ok]=check_io(model,graphics,-1,-1,[],[])
    if ok then
      rpar=[pitch;focus];
      ipar=[axes(:);finpix(:);length(fname);str2code(fname)]
      model.rpar=rpar; model.ipar=ipar;
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
case 'define' then
  fname='void'
  model=scicos_model()
  model.sim=list('hartmann',4)
  model.in=-1
  model.out=-1
  model.rpar=[0.000252;0.012246]
  model.ipar=[16;16;32;32;length(fname);str2code(fname)]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(axes);string(pitch);string(focus);string(finpix);fname]
  gr_i=['txt=['' Hartmann''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([3 2],model,exprs,gr_i)
end
endfunction

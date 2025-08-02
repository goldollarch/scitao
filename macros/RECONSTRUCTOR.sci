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

function [x,y,typ]=RECONSTRUCTOR(job,arg1,arg2)
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
    [ok,nsub,et,sit,fname,exprs]=getvalue(..
    ['Set Arroyo least squares reconstructors  parameter'],..
	    ['the number of subapertures';
	    'eigenvalue threshold';
	    'subaperture illumination threshold';
	    'reconstructor file name'],..
	    list('vec',1,'vec',1,'vec',1,'str',1),exprs)
   if ~ok then break,end
    
   fname=stripblanks(fname)
    
   [model,graphics,ok]=check_io(model,graphics,-1,[-1;-1],[],[])
    if ok then
      rpar=[et;sit];
      ipar=[nsub;length(fname);str2code(fname)];
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
  model.sim=list('reconstructor',4)
  model.in=-1
  model.out=[-1;-1]
  model.rpar=[1e-6;0.5]
  fname="void"
  model.ipar=[16;length(fname);str2code(fname)]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(nsub);string(et);string(sit);fname]
  gr_i=['txt=[''REC ''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([.5 2],model,exprs,gr_i)
end
endfunction

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

function [x,y,typ]=START(job,arg1,arg2)
// Copyright INRIA
x=[];y=[];typ=[];
select job
case 'plot' then //normal  position
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
    [ok,gd,gs,wl,exprs]=getvalue(..
      ['Set begin parameter'],..
	    ['the grid dimension';
	    'the grid size (meter)';
	    'wavelength (meter)'],..
	    list('vec',-1,'vec',1,'vec',1),exprs)
    if ~ok then break,end
    mess=[]
    if size(gd,'*')<>0 &size(gd,'*')<>2  then
      mess=[mess;'the grid dimension must be a 2 vector ';' ']
      ok=%f
    end
    if gs<0 then
      message('the grid size must be positive')
      ok=%f
    end
    if wl<0 then
      message('the wavelength must be positive')
      ok=%f
    end
    [model,graphics,ok]=check_io(model,graphics,[],1,1,[])
    if ok then
      rpar=[gs;wl];ipar=gd(:)
      model.rpar=rpar;
      model.ipar=ipar
      graphics.exprs=exprs;
      nout=2*gd(1)*gd(2)+30
      model.out=nout
      x.graphics=graphics;
      x.model=model
      break
    end
  end
  case 'define' then
  model=scicos_model()
  model.sim=list('start',4)
  model.in=[]
  model.out=-1
  model.evtin=1
  model.rpar=[0.01;1e-6]
  model.ipar=[256;256]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(gd);string(gs);string(wl)]
  gr_i=['txt=[''Begin ''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([3 2],model,exprs,gr_i)
end
endfunction

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

function [x,y,typ]=ZERNIKE(job,arg1,arg2)
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
    [ok,n,m,R,A,exprs]=getvalue(['Set the Zernike polynomials parameter'],..
	    ['the radial order n ';
	    'the azimuthal order m ';
	    'the radius R "
	    'the amplitude of aberration in R '],..
	    list('vect',1,'vect',1,'vec',1,'vec',1),exprs)
    if ~ok then break,end
    if n<0 then
      message('n must be positive')
      ok=%f
    end
    if abs(m)>n then
      message(' the absolution of m must be less or equal than n')
      ok=%f
    end
    if modulo(n,2)<>modulo(abs(m),2) then
      message(' n and m must be odd or even at the same time')
      ok=%f
    end
    [model,graphics,ok]=check_io(model,graphics,-1,-1,[],[])
    if ok then
      rpar=[R;A];ipar=[n;m]
      model.rpar=rpar;model.ipar=[n;m]
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
case 'define' then
  model=scicos_model()
  model.sim=list('zernike',4)
  model.in=-1
  model.out=-1
  model.rpar=[0.005;10]
  model.ipar=[2;2]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(n);
  string[m];
  string(R);
  string(A)]
  gr_i=['txt=['' Zernike ''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([3 2],model,exprs,gr_i)
end
endfunction

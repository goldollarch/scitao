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

function [x,y,typ]=DIAGNOSE(job,arg1,arg2)
// Copyright INRIA
x=[];y=[];typ=[]
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
  //dstate=model.in
  while %t do
    [ok,fname,heritance,exprs]=getvalue(..
	    'Set diagnose file parameters',..
	    ['file name';	'Accept herited events 0/1'],..
	    list('str',1,'vec',1),exprs)
    if ~ok then break,end //user cancel modification

    fname=stripblanks(fname)
    mess=[]
    if ~ok then
      message(['Some specified values are inconsistent:';
	         ' ';mess])
    end
    if ok then
      [model,graphics,ok]=check_io(model,graphics,-1,-1,ones(1-heritance,1),[])
    end
    
    if ok then
      ipar=[length(fname);str2code(fname)]
      model.dstate=[];
      
      model.ipar=ipar
      model.evtin=ones(1-heritance,1)
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
case 'define' then
  fname='outfile'
  model=scicos_model()
  model.sim=list('diagnose',4)
  model.in=-1
  model.out=-1
  model.evtin=1
  model.ipar=[length(fname);str2code(fname)]
  model.blocktype='c'
  model.dep_ut=[%t %f]

  exprs=[fname; string(0)'];
  gr_i=['txt=[''diagnose''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'')']
  x=standard_define([2 2],model,exprs,gr_i)
end
endfunction

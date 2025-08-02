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

function [x,y,typ]=ATMOSPHERE(job,arg1,arg2)
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
    [ok,nport,mtype,nlayers,pixsl,hvw,..
    r0,ggsmp,heighs,ps,subm,hwm,lfs,lwv,ppt,tm,fname,exprs]=getvalue(..
    ['Set refractive atmosphere model parameters'],..
	    ['the number of ports';
	    'atmosphere model type';
	    'layer number ';
	    'layer pixscales ';
	    'Hufnagel Valley pseudowind';
	    'Fred number r_0 and wavelength';
	    'Gemini GLAO study model parameters';
	    'layer heights ';
	    'layer power spectrum';
	    'subharmonic method parameters';
	    'Hardy wind model parameters';
	    'layer_foreshortening';
	    'layer_axes_wind_vector_aligned';
	    'propagation plan type';
	    'total simulation time';
	    'refractive layers file name'],..
	    list('vec',1,'vec',1,'vec',1,'vec',1,'vec',1,..
	    'vec',-1,'vec',-1,'vec',-1,'vec',-1,'vec',-1,'vec',-1,..
	    'vec',1,'vec',1,'vec',1,'vec',1,'str',1),exprs)
    if ~ok then break,end
    
    fname=stripblanks(fname)
    
    if nport==2 then
      [model,graphics,ok]=check_io(model,graphics,[-1;-1],[-1;-1],[],[])
    else
      [model,graphics,ok]=check_io(model,graphics,-1,-1,[],[])
   end
    
    if ok then
      rpar=[pixsl;hvw;r0(:);ggsmp(4);ps(2);ps(4);ps(5);hwm(:);tm;heighs(:)];
      ipar=[mtype;nlayers;ggsmp(1);ggsmp(2);ggsmp(3);ps(1);ps(3);
         subm(:);lfs;lwv;ppt;length(fname);str2code(fname)];
      model.rpar=rpar;model.ipar=ipar;
      graphics.exprs=exprs;
      x.graphics=graphics;
      x.model=model
      break
    end
  end
  
case 'define' then
  model=scicos_model()
  model.sim=list('atmosphere',4)
  model.in=-1;
  model.out=-1;
  
  hvw=21;
  r0=[0.05;0.5e-6];
  heighs=[1000]
  pixsl=0.02;
  hwm=[8;30;12000;4000]; 
  model.rpar=[pixsl;hvw;r0;10;-11/3.0;0.04;10;hwm;1;heighs]
  subm=[0;5;5;3]; 
  fname="void"
  model.ipar=[0;1;1;1;0;0;0;subm;0;0;0;length(fname);str2code(fname)]
  model.blocktype='c'
  model.dep_ut=[%f %t]
  
  exprs=[string(nport);string(mtype);
     string(nlayers); string(pixsl);
     string(hvw);string(r0);string(ggsmp);
     string(heighs);string(ps);string(subm);
     string(hwm);string(lfs);string(lwv);
     string(ppt);string(tm);fname]
  
  gr_i=['txt=['' ATM ''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'');']
  x=standard_define([.5 2],model,exprs,gr_i)
end
endfunction

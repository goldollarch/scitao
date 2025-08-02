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

function [x,y,typ]=DISPLAY(job,arg1,arg2)
// Copyright INRIA
x=[];y=[];typ=[]
select job
case 'plot' then
  standard_draw(arg1)
case 'getinputs' then
  [x,y,typ]=standard_inputs(o)
case 'getoutputs' then
  x=[];y=[];typ=[];
case 'getorigin' then
  [x,y]=standard_origin(arg1)
case 'set' then
  x=arg1;
  graphics=arg1.graphics;exprs=graphics.exprs
  if size(exprs)<5 then exprs(5)=emptystr(),end // compatibility
  model=arg1.model;
  //dstate=model.in
  while %t do
    [ok,typ,win,wpos,wdim,heritance,nom,exprs]=getvalue(..
	'Set Scope parameters',..
	['display type';
       'window ID (-1 for automatic)';
	'Output window position';
	'Output window sizes';
	'Accept herited events (0/1)';
	'Name of Scope (label&Id)'],..
	 list('vec',1,'vec',1,'vec',-1,'vec',-1,'vec',1,'str',1),..
	 exprs)
    if ~ok then break,end //user cancel modification
    mess=[]
    if size(wpos,'*')<>0 &size(wpos,'*')<>2 then
      mess=[mess;'Window position must be [] or a 2 vector';' ']
      ok=%f
    end
    if size(wdim,'*')<>0 &size(wdim,'*')<>2 then
      mess=[mess;'Window dim must be [] or a 2 vector';' ']
      ok=%f
    end
    if win<-1 then
      mess=[mess;'Window number can''t be  < -1';' ']
      ok=%f
    end
    if ~or(heritance==[0 1]) then
      mess=[mess;'Accept herited events must be 0 or 1';' ']
      ok=%f
    end
    if ~ok then
      message(['Some specified values are inconsistent:';
	         ' ';mess])
	   end
    if ok then
      [model,graphics,ok]=check_io(model,graphics,-1,[],ones(1-heritance,1),[])
    end
    
    if ok then
      if wpos==[] then wpos=[-1;-1];end
      if wdim==[] then wdim=[-1;-1];end
      ipar=[typ;win;wpos(:);wdim(:)]
      model.dstate=[];
      
      model.ipar=ipar
      model.evtin=ones(1-heritance,1)
      model.label=nom;
      graphics.id=nom
      graphics.exprs=exprs;
      x.graphics=graphics;x.model=model
      break
    end
  end
case 'define' then
  win=-1;
  wdim=[400;400]
  wpos=[-1;-1]

  model=scicos_model()
  model.sim=list('display',4)
  model.in=-1
  model.evtin=1
  model.ipar=[0;win;wpos;wdim]
  model.blocktype='c'
  model.dep_ut=[%t %f]

  exprs=[string(typ);
       string(win);
	sci2exp([]);
	sci2exp(wdim);
       string(0)';
	emptystr()]; //label-id
  gr_i=['txt=[''display ''];';
    'xstringb(orig(1),orig(2),txt,sz(1),sz(2),''fill'')']
  x=standard_define([2 2],model,exprs,gr_i)
end
endfunction


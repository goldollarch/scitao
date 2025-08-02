
mode(-1);
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

m=1;
nm=1e-9*m;
mm=1e-3*m;
cm=1e-2*m;

lambda=550*nm;
sizes=5*mm;
N=100;

R=0.12*mm;
d=0.5*mm;
z=5*cm;

F=begin(sizes,lambda,N);
F1=circ_ap(F,R,d,0);
F2=circ_ap(F,R,-d,0);
F=b_mix(F1,F2);

xbasc();
drawlater();

ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);

for i=1:16
   F=forvard(F,z);
   E=field_int(F);
   subplot(4,4,i);
   grayplot(1:N,1:N,E,strf="030",rect=[1,1,N,N]);
   Str=sprintf('z=%4.1f cm',i*z/cm)
   xtitle(Str);
end

drawnow();


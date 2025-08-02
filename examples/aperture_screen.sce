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

lambda=1000*nm;
sizes=10*mm;
N=200;

R=2.5*mm;

x=1:N;y=1:N;
F=begin(sizes,lambda,N);
H=rect_screen(F,0.001,0.001, -0.0015, -0.002, 0.);
G=rect_screen(H,0.001,0.0035, -0.002, 0.0025, 45);
I=field_int(G);

xbasc();
drawlater();

subplot(2,2,1);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('rectangular screen')
xtitle(Str);

Y=forvard(G,1);
I=field_int(Y);

subplot(2,2,2);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('screen diffraction')
xtitle(Str);

T=circ_screen(G,0.0005,0.0025,-0.003);
Y=circ_screen(T,0.0007,0.001,0.0015);
I=field_int(Y);

subplot(2,2,3);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('circular screen')
xtitle(Str);

Y=forvard(Y,1);
I=field_int(Y);
subplot(2,2,4);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);
Str=sprintf('screen difractive ')
xtitle(Str);

drawnow();

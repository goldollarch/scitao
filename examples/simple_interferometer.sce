
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

xbasc();
drawlater();
ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);

// two round pinholes
N=200;x=1:N;y=1:N;
F=begin(0.005,0.55e-6,N);

E1=circ_ap(F,0.00012,0.,-0.0005);
E2=circ_ap(F,0.00012,0.,0.0005);
Y=b_mix(E1,E2);
I=field_int(Y);
subplot(3,3,1);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,0.5);
I=field_int(G);
subplot(3,3,2);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,1);
I=field_int(G);
subplot(3,3,3);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

// three round holes
E1=circ_ap(F,0.00012,-0.00025,-0.0005);
E2=circ_ap(F,0.00012,-0.00025,0.0005);
E3=circ_ap(F,0.00012,0.00025,0.);
Y1=b_mix(E1,E2);
Y=b_mix(E3,Y1);
I=field_int(Y);
subplot(3,3,4);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,0.5);
I=field_int(G);
subplot(3,3,5);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,1);
I=field_int(G);
subplot(3,3,6);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

// two slits
E1=rect_ap(F,0.0025, 0.0001, 0., -0.0005, 0.);
E2=rect_ap(F,0.0025, 0.0001, 0., 0.0005, -15);
Y=b_mix(E1,E2);
I=field_int(Y);
subplot(3,3,7);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,0.5);
I=field_int(G);
subplot(3,3,8);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

G=forvard(Y,1);
I=field_int(G);
subplot(3,3,9);
grayplot(x,y,I,strf="030",rect=[1,1,N,N]);

drawnow();

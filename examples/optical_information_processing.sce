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

N=128;
x=1:N;y=1:N;

// 	create the intensity mask:
wavelength=1e-6;
pixscl=0.1;sz=N*pixscl;
dwf=begin(2*sz,wavelength,N);
ap=til_hex_ap(sz,sz/8,sz/13,sz/100);
dwf=aperture_transform(dwf,ap);
file_pgm(dwf,"test",-1,1);

// 	importing the intensity mask:
f1=begin(1e-2,5e-7,N);
f2=fil_ter(f1,"int","subs","test.pgm");
int1=field_int(f2);

f3=forvard(f2,1);

//filtered with a high-frequency component cut (center), 
//	lens, focal filter, free space...
f4=lens(f3,1);
f5=forvard(f4,1);

f6=circ_ap(f5,0.0015);
f7=forvard(f6,1);
f8=lens(f7,1);
f9=forvard(f8,1);
int2=field_int(f9);

//filtered with a low frequency components cut (center), 
//	lens, focal filter, free space...
f6=circ_screen(f5,0.0015);
f7=forvard(f6,1);
f8=lens(f7,1);
f9=forvard(f8,1);
int3=field_int(f9);

//filtered with other components cut , 
//	lens, focal filter, free space...
f6=gauss(f5,5e-4);
f7=forvard(f6,1);
f8=lens(f7,1);
f9=forvard(f8,1);
int4=field_int(f9);

subplot(2,2,1);
grayplot(x,y,int1,strf="030",rect=[1,1,N,N]);

subplot(2,2,2);
grayplot(x,y,int2,strf="030",rect=[1,1,N,N]);
Str=sprintf('filter high-frequency')
xtitle(Str);

subplot(2,2,3);
grayplot(x,y,int3,strf="030",rect=[1,1,N,N]);
Str=sprintf('filter low-frequency')
xtitle(Str);

subplot(2,2,4);
grayplot(x,y,int4,strf="030",rect=[1,1,N,N]);
Str=sprintf('filter other-frequency')
xtitle(Str);

drawnow();


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
// Simulation of a Michelson interferometer

m=1; 
nm=1e-9*m; 
mm=1e-3*m; 
cm=1e-2*m;

rad=1; 
mrad=1e-3*rad;

N=250;
sz=25*mm;
lambda=500*nm;

R=12*mm;
z1=10*cm;
z4=10*cm;
RBS=0.5;

ty=0.0*mrad;
tx=0.0*mrad;
//tx=0.2*mrad;

f=600*cm;
z3=30*cm;
z20=50*cm;

// A weak converging beam using a weak positive lens:
F=begin(sz,lambda,N);
F=circ_ap(F,R);
F=lens(F,f);

// Propagation to the beamsplitter:
F=forvard(F,z1);

// Splitting the beam and propagation to mirror #2:
F2=absorber(F,1-RBS);
F2=forvard(F2,z3);

// Introducing tilt and propagating back to the beamsplitter:
F2=tilt(F2,tx,ty);
F2=forvard(F2,z3);
F2=absorber(F2,RBS);

// Splitting off the second beam:
F10=absorber(F,RBS);

// Initializing the screen and the array for storage of the movie:

graph=scf(); 
clf(graph)
    
graph.pixmap='on'
Axe=graph.children
Axe.isoview='on'
    
x = 1:N;y = 1:N;
z = zeros(N,N);
  
ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);
colorbar(0,1)

Sgrayplot(x,y,z,strf="042",zminmax=[0,1])
    
Axe=graph.children
Axe.box='off'
    
c=gce();e=c.children

num=10;
// Scanning mirror #1:
for FRAME=1:num
  
  z2=z20+(lambda/num)*(FRAME-1);
  F1=forvard(F10,z2*2);
  F1=absorber(F1,1-RBS);
  
  // Recombining the two beams and propagation to the screen:
  F=b_mix(F1,F2);
  F=forvard(F,z4);
  
  I=field_int(F);
  
  e.data(:,3)=matrix(I,-1,1)
  show_pixmap();
  
end;


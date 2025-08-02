mode(-1)
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

//Twyman-Green shearing interferometer.

F=begin(0.03,5e-7);
H=circ_ap(F,0.005);
G=forvard(H,0.5);
[S1,S2]=b_split(G,0.3);

K1=forvard(S1,1);
E1=b_split(K1,0.7);

K2=forvard(S2,0.4);
T2=zernike(K2,3, 1, 0.005, 25);
Y2=forvard(T2,0.4);
E2=b_split(Y2,0.3);

E=b_mix(E1,E2);
K=forvard(E,1.);
S=interpol(K,0.012,128,0.,0.,0.,1.);

field_plot(S);
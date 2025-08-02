
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

lambda=1000*nm;
sizes=30*mm;
R=15*mm;

N=128;

F=begin(sizes,lambda,N);
header=field_contents(F);

intensity=zeros(N,N);
phase=zeros(N,N);

for i=1:N
	for j=1:N
		intensity(i,j)=abs(sin(i/10)*cos(j/5));
		phase(i,j)=cos(i/10)*sin(j/5);
	end
end

G=create_field(header,intensity,phase);
H=circ_ap(G,R);
field_plot(H);

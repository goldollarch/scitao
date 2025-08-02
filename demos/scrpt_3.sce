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

N=100;
f1=begin(0.01,1e-6,N);
f2=circ_ap(f1,0.005);

xbasc();
drawlater();

ncolor=60;
h=hotcolormap(ncolor);
xset('colormap',h);

n=0;
for i=0:6
	ii=modulo(i,2)
	for j=ii:2:i
		n=n+1;
		f3=zernike(f2,i,j,0.005,1);
		pha=field_pha(f3);
		subplot(4,4,n); 
		grayplot(1:N,1:N,pha,strf="030",rect=[1,1,N,N]);
    	Str=sprintf('n=%d, m=%d',i,j)
    	xtitle(Str);		
	end
end

drawnow();

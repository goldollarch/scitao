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

function scao_ui()

  ok=%f;
  fin=%t;

  scao_para=[]
 
  list_heights = []
  heights=strcat(list_heights,'|');
  
  list_objects =['2.2e-6' '8' '8']
  objects=strcat(list_objects,'|');
  
  f = figure("Position",[150 50 1000 820],...
	     "BackgroundColor",[0.7 0.9 0.4],...
	     "Unit", "pixel");

  m=uimenu(f,'label', 'menu');
  mhelp=uimenu(f,'label', 'help', 'callback', "scao_help()");
  mexit=uimenu(f,'label', 'exit', 'callback', "exit");

  txt0 = uicontrol(f, "Position",[90 760 100 40],...
	"Style","text","String","SCAO ","fontname",...
	"Times Bold Italic","fontsize",40,...
	"BackgroundColor",[0.7 0.9 0.4]);

  frame01 = uicontrol(f, "Position",[265 675 180 120],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt01 = uicontrol(f, "Position",[280 765 100 30],...
	"Style","text","String","time domain","fontsize",18,...
	"BackgroundColor",[0.9 0.9 0.9]);
  txt011 = uicontrol(f, "Position",[280 735 70 20],...
	"Style","text","String","nsteps","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit011 = uicontrol(f, "Position" ,[360 735 60 20],...
	"Style","edit","String","25",...
	"BackgroundColor",[1 1 1]);		
  txt012 = uicontrol(f, "Position",[280 700 70 20],...
	"Style","text","String","timestep","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit012 = uicontrol(f, "Position" ,[360 700 60 20],...
	"Style","edit","String","0.002",...
	"BackgroundColor",[1 1 1]);		
	
  frame02 = uicontrol(f, "Position",[465 675 180 120],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt02 = uicontrol(f, "Position",[475 765 80 30],...
	"Style","text","String","other time",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt021 = uicontrol(f, "Position",[470 747 100 20],...
	"Style","text","String","prop delay","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit021 = uicontrol(f, "Position" ,[570 750 60 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt022 = uicontrol(f, "Position",[470 720 100 20],...
	"Style","text","String","prop interval","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit022 = uicontrol(f, "Position" ,[570 720 60 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt023 = uicontrol(f, "Position",[470 690 100 20],...
	"Style","text","String","dm loop delay","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit023 = uicontrol(f, "Position" ,[570 690 60 20],...
	"Style","edit","String","5",...
	"BackgroundColor",[1 1 1]);		
	
  frame03 = uicontrol(f, "Position",[660 675 140 120],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt03 = uicontrol(f, "Position",[670 765 100 30],...
	"Style","text","String","propagation" ,...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt031 = uicontrol(f, "Position",[675 735 70 20],...
	"Style","text","String","geometry","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit031 = uicontrol(f, "Position",[750 735 30 20],...
	"Style", "edit","String","1",...
	"BackgroundColor",[1 1 1] );
  txt032 = uicontrol(f, "Position",[675 700 70 20],...
	"Style","text","String","diffractive","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit032 = uicontrol(f, "Position",[750 700 30 20],...
	"Style", "edit","String","0",...
	"BackgroundColor",[1 1 1] );
		
  frame04 = uicontrol(f, "Position",[820 675 160 120],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt03 = uicontrol(f, "Position",[830 765 120 30],...
	"Style","text","String","output verbose" ,...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt041 = uicontrol(f, "Position",[835 735 70 20],...
	"Style","text","String","verbose","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit041 = uicontrol(f, "Position" ,[915 735 40 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt042 = uicontrol(f, "Position",[835 700 70 20],...
	"Style","text","String","vverbose","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit042 = uicontrol(f, "Position" ,[915 700 40 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
	
  frame10= uicontrol(f, "Position",[20 675 225 70],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt10 = uicontrol(f, "Position",[35 720 60 25],...
	"Style","text","String","beacon",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt101 = uicontrol(f, "Position",[50 695 70 20],...
	"Style","text","String","wavelengh","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit101 = uicontrol(f, "Position" ,[130 695 80 25],...
	"Style","edit","String","0.6e-6",...
	"BackgroundColor",[1 1 1]);		
	
  frame11= uicontrol(f, "Position",[20 350 225 300],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt11 = uicontrol(f, "Position",[35 625 60 25],...
	"Style","text","String","objects",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt110 = uicontrol(f, "Position",[30 590 40 20],...
	"Style","text","String","num","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit110 = uicontrol(f, "Position" ,[70 590 40 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt111 = uicontrol(f, "Position",[30 555 40 20],...
	"Style","text","String","nx","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit111 = uicontrol(f, "Position" ,[70 555 40 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt112 = uicontrol(f, "Position",[30 520 40 20],...
	"Style","text","String","ny","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit112 = uicontrol(f, "Position" ,[70 520 40 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt113 = uicontrol(f, "Position",[30 485 40 20],...
	"Style","text","String","gap","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit113 = uicontrol(f, "Position" ,[70 485 40 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
	
  frame112= uicontrol(f, "Position",[30 360 200 90],...
	"Style","frame","BackgroundColor",[0.7 0.7 0.7]);
  txt114 = uicontrol(f, "Position",[40 425 50 20],...
	"Style","text","String","lamda",...
	"fontsize",15,"BackgroundColor",[0.7 0.7 0.7]);
  edit114 = uicontrol(f, "Position" ,[40 400 50 20],...
	"Style","edit","String","2.2e-6",...
	"BackgroundColor",[1 1 1]);		
  txt115 = uicontrol(f, "Position",[105 425 50 20],...
	"Style","text","String","img_sz",...
	"fontsize",15,"BackgroundColor",[0.7 0.7 0.7]);
  edit115 = uicontrol(f, "Position" ,[105 400 50 20],...
	"Style","edit","String","8",...
	"BackgroundColor",[1 1 1]);		
  txt116 = uicontrol(f, "Position",[170 425 50 20],...
	"Style","text","String","sample",...
	"fontsize",15,"BackgroundColor",[0.7 0.7 0.7]);
  edit116 = uicontrol(f, "Position" ,[170 400 50 20],...
	"Style","edit","String","8",...
	"BackgroundColor",[1 1 1]);		
  list_obj = uicontrol(f, "Position" ,[130 470 90 150],...
	"Style", "listbox", "String" ,objects,...
	"BackgroundColor",[1 1 1]);
  add11 = uicontrol(f, "Position",[70 370 50 20],...
	"Style","pushbutton","String","add",...
	"callback"  , "add_list_objects()" ,...
	"BackgroundColor",[0 1 0]);  	
  del11 = uicontrol(f, "Position",[140 370 50 20],...
	"Style","pushbutton","String","delete",...
	"callback"  , "del_list_objects()" ,...
	"BackgroundColor",[0 1 0]);  	

  frame12= uicontrol(f, "Position",[265 350 225 300],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt12 = uicontrol(f, "Position",[280 620 80 30],...
	"Style","text","String","aperture","fontsize",18,...
	"BackgroundColor",[0.9 0.9 0.9]);
  txt121 = uicontrol(f, "Position",[295 595 70 20],...
	"Style","text","String","type","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit121 = uicontrol(f, "Position" ,[385 595 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt122 = uicontrol(f, "Position",[280 560 100 20],...
	"Style","text","String","inner diameter","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit122 = uicontrol(f, "Position" ,[385 560 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt123 = uicontrol(f, "Position",[280 525 100 20],...
	"Style","text","String","outer diameter","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit123 = uicontrol(f, "Position" ,[385 525 80 20],...
	"Style","edit","String","2",...
	"BackgroundColor",[1 1 1]);		
  txt124 = uicontrol(f, "Position",[280 490 100 20],...
	"Style","text","String","nspiders","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit124 = uicontrol(f, "Position" ,[385 490 80 20],...
	"Style","edit","String","4",...
	"BackgroundColor",[1 1 1]);		
  txt125 = uicontrol(f, "Position",[280 455 100 20],...
	"Style","text","String","spider width","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit125 = uicontrol(f, "Position" ,[385 455 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt126 = uicontrol(f, "Position",[280 420 100 20],...
	"Style","text","String","edge length","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit126 = uicontrol(f, "Position" ,[385 420 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt127 = uicontrol(f, "Position",[280 385 100 20],...
	"Style","text","String","gap size","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit127 = uicontrol(f, "Position" ,[385 385 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
	
  frame21= uicontrol(f, "Position",[510 510 225 140],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt21 = uicontrol(f, "Position",[525 620 130 30],...
	"Style","text","String","TTM/controllers",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt211 = uicontrol(f, "Position",[525 595 100 20],...
	"Style","text","String","velocity","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit211 = uicontrol(f, "Position" ,[630 595 80 20],...
	"Style","edit","String","100",...
	"BackgroundColor",[1 1 1]);		
  txt212 = uicontrol(f, "Position",[525 560 100 20],...
	"Style","text","String","proportional","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit212 = uicontrol(f, "Position" ,[630 560 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt213 = uicontrol(f, "Position",[525 525 100 20],...
	"Style","text","String","integral","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit213 = uicontrol(f, "Position" ,[630 525 80 20],...
	"Style","edit","String","5e-2",...
	"BackgroundColor",[1 1 1]);		

  frame22= uicontrol(f, "Position",[510 350 225 140],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt22 = uicontrol(f, "Position",[525 460 130 30],...
	"Style","text","String","DM/controllers",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt221 = uicontrol(f, "Position",[525 435 100 20],...
	"Style","text","String","velocity","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit221 = uicontrol(f, "Position" ,[630 435 80 20],...
	"Style","edit","String","2.67e-4",...
	"BackgroundColor",[1 1 1]);		
  txt222 = uicontrol(f, "Position",[525 400 100 20],...
	"Style","text","String","proportional","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit222 = uicontrol(f, "Position" ,[630 400 80 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt223 = uicontrol(f, "Position",[525 365 100 20],...
	"Style","text","String","integral","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit223 = uicontrol(f, "Position" ,[630 365 80 20],...
	"Style","edit","String","1e-7",...
	"BackgroundColor",[1 1 1]);		

  frame23= uicontrol(f, "Position",[755 470 225 180],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt23 = uicontrol(f, "Position",[770 620 150 30],...
	"Style","text","String","Hartmann sensor",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  txt231 = uicontrol(f, "Position",[770 595 100 20],...
	"Style","text","String","focal length","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit231 = uicontrol(f, "Position" ,[875 595 80 20],...
	"Style","edit","String","0.012246",...
	"BackgroundColor",[1 1 1]);		
  txt1232 = uicontrol(f, "Position",[770 560 100 20],...
	"Style","text","String","lenslet pitch","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit232 = uicontrol(f, "Position" ,[875 560 80 20],...
	"Style","edit","String","0.000252",...
	"BackgroundColor",[1 1 1]);		
  txt233 = uicontrol(f, "Position",[770 525 100 20],...
	"Style","text","String","pixels / lenslet","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit233 = uicontrol(f, "Position" ,[875 525 80 20],...
	"Style","edit","String","32",...
	"BackgroundColor",[1 1 1]);		
  txt234 = uicontrol(f, "Position",[770 490 100 20],...
	"Style","text","String","pixels / xform","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit234 = uicontrol(f, "Position" ,[875 490 80 20],...
	"Style","edit","String","32",...
	"BackgroundColor",[1 1 1]);		

  frame3= uicontrol(f, "Position",[20 30 715 300],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt3 = uicontrol(f, "Position",[35 300 160 30],...
	"Style","text","String","atmospheric models",...
	"fontsize",18,"BackgroundColor",[0.9 0.9 0.9]);
  frame31= uicontrol(f, "Position",[35 45 330 120],...
	"Style","frame","BackgroundColor",[0.7 0.7 0.7]);
  txt31 = uicontrol(f, "Position",[50 145 150 20],...
	"Style","text","String","subharmonic method",...
	"fontsize",16,"BackgroundColor",[0.7 0.7 0.7]);
  txt311 = uicontrol(f, "Position",[60 110 50 20],...
	"Style","text","String","type","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit311 = uicontrol(f, "Position" ,[120 110 50 20],...
	"Style","edit","String","2",...
	"BackgroundColor",[1 1 1]);		
  txt312 = uicontrol(f, "Position",[60 70 50 20],...
	"Style","text","String","depth","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit312 = uicontrol(f, "Position" ,[120 70 50 20],...
	"Style","edit","String","3",...
	"BackgroundColor",[1 1 1]);		
  txt313 = uicontrol(f, "Position",[180 110 100 20],...
	"Style","text","String","subpixels/level","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit313 = uicontrol(f, "Position" ,[290 110 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt314 = uicontrol(f, "Position",[180 70 100 20],...
	"Style","text","String","subpixels/pixel","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit314 = uicontrol(f, "Position" ,[290 70 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  frame32= uicontrol(f, "Position",[35 175 330 120],...
	"Style","frame","BackgroundColor",[0.7 0.7 0.7]);
  txt32 = uicontrol(f, "Position",[50 275 130 20],...
	"Style","text","String","Hardy wind model",...
	"fontsize",16,"BackgroundColor",[0.7 0.7 0.7]);
  txt321 = uicontrol(f, "Position",[60 240 50 20],...
	"Style","text","String","height","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit321 = uicontrol(f, "Position" ,[120 240 50 20],...
	"Style","edit","String","10000",...
	"BackgroundColor",[1 1 1]);		
  txt322 = uicontrol(f, "Position",[50 200 70 20],...
	"Style","text","String","thickness","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit322 = uicontrol(f, "Position" ,[120 200 50 20],...
	"Style","edit","String","5000",...
	"BackgroundColor",[1 1 1]);		
  txt323 = uicontrol(f, "Position",[180 240 100 20],...
	"Style","text","String","ground wind","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit323 = uicontrol(f, "Position" ,[290 240 50 20],...
	"Style","edit","String","5",...
	"BackgroundColor",[1 1 1]);		
  txt324 = uicontrol(f, "Position",[180 200 100 20],...
	"Style","text","String","troposph wind","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit324 = uicontrol(f, "Position" ,[290 200 50 20],...
	"Style","edit","String","30",...
	"BackgroundColor",[1 1 1]);		
  frame33= uicontrol(f, "Position",[380 45 340 270],...
	"Style","frame","BackgroundColor",[0.7 0.7 0.7]);
  txt33 = uicontrol(f, "Position",[400 295 210 20],...
	"Style","text","String","refractive atmospheric model",...
	"fontsize",16,"BackgroundColor",[0.7 0.7 0.7]);
  txt330 = uicontrol(f, "Position",[420 270 50 20],...
	"Style","text","String","type","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit330 = uicontrol(f, "Position" ,[480 270 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt331 = uicontrol(f, "Position",[570 270 50 20],...
	"Style","text","String","nlayers","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit331 = uicontrol(f, "Position" ,[630 270 50 20],...
	"Style","edit","String","0",...
	"BackgroundColor",[1 1 1]);		
  txt332 = uicontrol(f, "Position",[400 240 120 20],...
	"Style","text","String","r_0_meters","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit332 = uicontrol(f, "Position" ,[530 240 50 20],...
	"Style","edit","String","0.05",...
	"BackgroundColor",[1 1 1]);		
  txt333 = uicontrol(f, "Position",[400 210 120 20],...
	"Style","text","String","wavelength","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit333 = uicontrol(f, "Position" ,[530 210 50 20],...
	"Style","edit","String","0.5e-6",...
	"BackgroundColor",[1 1 1]);		
  txt334 = uicontrol(f, "Position",[400 180 120 20],...
	"Style","text","String","pseudowind","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit334 = uicontrol(f, "Position" ,[530 180 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt335 = uicontrol(f, "Position",[400 150 120 20],...
	"Style","text","String","outer scale","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit335 = uicontrol(f, "Position" ,[530 150 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt336 = uicontrol(f, "Position",[400 120 120 20],...
	"Style","text","String","ground quality","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit336 = uicontrol(f, "Position" ,[530 120 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt337 = uicontrol(f, "Position",[400 90 120 20],...
	"Style","text","String","anisoplan quality","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit337 = uicontrol(f, "Position" ,[530 90 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt338 = uicontrol(f, "Position",[400 60 120 20],...
	"Style","text","String","extended profile","fontsize",15,...
	"BackgroundColor",[0.7 0.7 0.7]);
  edit338 = uicontrol(f, "Position" ,[530 60 50 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt339 = uicontrol(f, "Position",[630 240 50 20],...
	"Style","text","String","heights",...
	"fontsize",15,"BackgroundColor",[0.7 0.7 0.7]);
  list_hgh = uicontrol(f, "Position" ,[610 130 90 110],...
	"Style", "listbox", "String" ,heights,...
	"BackgroundColor",[1 1 1]);
  edit339 = uicontrol(f, "Position" ,[630 100 50 20],...
	"Style","edit","String","",...
	"BackgroundColor",[1 1 1]);		
  add33 = uicontrol(f, "Position",[600 65 50 20],...
	"Style","pushbutton","String","add",...
	"callback"  , "add_list_heights()" ,...
	"BackgroundColor",[0 1 0]);  	
  del33 = uicontrol(f, "Position",[660 65 50 20],...
	"Style","pushbutton","String","delete",...
	"callback"  , "del_list_heights()" ,...
	"BackgroundColor",[0 1 0]);  	

  frame4= uicontrol(f, "Position",[755 270 225 180],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt4 = uicontrol(f, "Position",[770 420 110 30],...
	"Style","text","String","reconstructor","fontsize",18,...
	"BackgroundColor",[0.9 0.9 0.9]);
  txt41 = uicontrol(f, "Position",[770 395 100 20],...
	"Style","text","String","nsubaps","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit41 = uicontrol(f, "Position" ,[875 395 80 20],...
	"Style","edit","String","16",...
	"BackgroundColor",[1 1 1]);		
  txt42 = uicontrol(f, "Position",[770 360 100 20],...
	"Style","text","String","real threshold","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit42 = uicontrol(f, "Position" ,[875 360 80 20],...
	"Style","edit","String","0.5",...
	"BackgroundColor",[1 1 1]);		
  txt43 = uicontrol(f, "Position",[770 325 100 20],...
	"Style","text","String","eigen threshold","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit43 = uicontrol(f, "Position" ,[875 325 80 20],...
	"Style","edit","String","1e-6",...
	"BackgroundColor",[1 1 1]);		
  txt44 = uicontrol(f, "Position",[770 290 100 20],...
	"Style","text","String","filename","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit44 = uicontrol(f, "Position" ,[875 290 80 20],...
	"Style","edit","String","void",...
	"BackgroundColor",[1 1 1]);		

  frame5= uicontrol(f, "Position",[755 100 225 150],...
	"Style","frame","BackgroundColor",[0.9 0.9 0.9]);
  txt5 = uicontrol(f, "Position",[770 220 50 30],...
	"Style","text","String","others","fontsize",18,...
	"BackgroundColor",[0.9 0.9 0.9]);
  txt51 = uicontrol(f, "Position",[770 195 120 20],...
	"Style","text","String","random seed","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit51 = uicontrol(f, "Position" ,[895 195 60 20],...
	"Style","edit","String","1",...
	"BackgroundColor",[1 1 1]);		
  txt52 = uicontrol(f, "Position",[770 160 120 20],...
	"Style","text","String","wavefront scale","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit52 = uicontrol(f, "Position" ,[895 160 60 20],...
	"Style","edit","String","0.02",...
	"BackgroundColor",[1 1 1]);		
  txt53 = uicontrol(f, "Position",[770 125 120 20],...
	"Style","text","String","layer scale","fontsize",15,...
	"BackgroundColor",[0.9 0.9 0.9]);
  edit53 = uicontrol(f, "Position" ,[895 125 60 20],...
	"Style","edit","String","0.02",...
	"BackgroundColor",[1 1 1]);		

  bcheck = uicontrol(f, "Position" ,[780 40 60 30],...
	"Style", "pushbutton","String", "check",...
	"callback"  , "check_and_init()","BackgroundColor",[0 1 0] );
  bexecute = uicontrol(f, "Position",[850 40 50 30],...
	"Style", "pushbutton","String", "execute",...
	"callback" ,"do_execute()","BackgroundColor",[1 0 0] );  		 
  bquit = uicontrol(f, "Position",[910 40 50 30],...
	"Style", "pushbutton","String", "cancel",...
	"callback" ,"close(f)","BackgroundColor",[0 1 0] );  		 

  while fin,end
  
//  fin=%f
  while ~fin
    sleep(1)
    if findobj('label','menu')==[] then
      return;
    end
  end
  
  close(f)
  return
  
endfunction

//////////////////////////////

function scao_help()
	help scao
endfunction

function add_list_heights()
  if edit339<>0&list_hgh<>0 then
    new = get(edit339,'String'); 
    if (list_heights == '') then
      list_heights=[];
    end
    list_heights = [list_heights  new];
    set(list_hgh,'String',strcat(list_heights,'|'));
    list_heights=resume(list_heights)
  end
endfunction

function del_list_heights()
  if edit339<>0&list_hgh<>0 then
    del = get(list_hgh,'Value');
    if (list_heights == '') then
      list_heights=[];
    end
    list_heights =[list_heights(1:del-1) list_heights(del+1:$)];
    set(list_hgh,'String',strcat(list_heights,'|'));
    list_heights=resume(list_heights);
  end
endfunction

function add_list_objects()
  if edit114<>0&edit115<>0&edit116<>0&list_obj<>0 then
    new1 = get(edit114,'String'); 
    new2 = get(edit115,'String'); 
    new3 = get(edit116,'String'); 
    if (list_objects == '') then
      list_objects=[];
    end
    list_objects = [list_objects  new1 new2 new3];
    set(list_obj,'String',strcat(list_objects,'|'));
    list_objects=resume(list_objects)
  end
endfunction

function del_list_objects()
  if edit114<>0&list_obj<>0 then
    del = get(list_obj,'Value');
    if (list_objects == '') then
      list_objects=[];
    end
    list_objects =[list_objects(1:del-1) list_objects(del+3:$)];
    set(list_obj,'String',strcat(list_objects,'|'));
    list_objects=resume(list_objects);
  end
endfunction

///////////////////////////////////////////////////////

function check_and_init()
  
  fin=%f;  
  
  nstep=get(edit011,'String');
  if (nstep == '') then
      x_message(['the parameter nsteps is not valid']); 
      fin=%t;
  else
     execstr('nsteps = '+nstep);
     scao_para=[scao_para nsteps];
  end
  
  dm_loop_delay=get(edit023 ,'String');
  if (dm_loop_delay == '') then
      x_message(['the parameter nsteps is not valid']); 
      fin=%t;
  else
     execstr('dm_loop_delay = '+dm_loop_delay);
     scao_para=[scao_para dm_loop_delay];
  end
  
  prop_interval=get(edit022 ,'String');
  if (prop_interval == '') then
      x_message(['the parameter nsteps is not valid']); 
      fin=%t;
  else
     execstr('prop_interval = '+prop_interval);
     scao_para=[scao_para prop_interval];
  end
  
  prop_delay=get(edit021 ,'String');
  if (prop_delay == '') then
      x_message(['the parameter nsteps is not valid']); 
      fin=%t;
  else
     execstr('prop_delay = '+prop_delay);
     scao_para=[scao_para prop_delay];
  end
  
  timestep=get(edit012 ,'String');
  if (timestep == '') then
      x_message(['the parameter nsteps is not valid']); 
      fin=%t;
  else
     execstr('timestep = '+timestep);
     scao_para=[scao_para timestep];
  end

  rand_seed=get(edit51,'String');
  if (rand_seed == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('rand_seed = '+rand_seed);
     scao_para=[scao_para rand_seed];
  end

  ttm_vel=get(edit211,'String');
  if (ttm_vel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('ttm_vel = '+ttm_vel);
     scao_para=[scao_para ttm_vel];
  end

  ttm_prop=get(edit212,'String');
  if (ttm_prop == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('ttm_prop = '+ttm_prop);
     scao_para=[scao_para ttm_prop];
  end

  ttm_int=get(edit213,'String');
  if (ttm_int == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('ttm_int = '+ttm_int);
     scao_para=[scao_para ttm_int];
  end

  dm_vel=get(edit221,'String');
  if (dm_vel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('dm_vel = '+dm_vel);
     scao_para=[scao_para dm_vel];
  end

  dm_prop=get(edit222,'String');
  if (dm_prop == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('dm_prop = '+dm_prop);
     scao_para=[scao_para dm_prop];
  end

  dm_int=get(edit223,'String');
  if (dm_int == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('dm_int = '+dm_int);
     scao_para=[scao_para dm_int];
  end

  beacon_wvl=get(edit101,'String');
  if (beacon_wvl == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('beacon_wvl = '+beacon_wvl);
     scao_para=[scao_para beacon_wvl];
  end

  objects_num=get(edit110,'String');
  if (objects_num == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('objects_num = '+objects_num);
     scao_para=[scao_para objects_num];
  end

  subm_type=get(edit311,'String');
  if (subm_type == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('subm_type = '+subm_type);
     scao_para=[scao_para subm_type];
  end

  subm_depth=get(edit312,'String');
  if (subm_depth == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('subm_depth = '+subm_depth);
     scao_para=[scao_para subm_depth];
  end

  subpix_level=get(edit313,'String');
  if (subpix_level == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('subpix_level = '+subpix_level);
     scao_para=[scao_para subpix_level];
  end

  subpix_pix=get(edit314,'String');
  if (subpix_pix == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('subpix_pix = '+subpix_pix);
     scao_para=[scao_para subpix_pix];
  end

  wind_high=get(edit321,'String');
  if (wind_high == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('wind_high = '+wind_high);
     scao_para=[scao_para wind_high];
  end

  wind_thick=get(edit322,'String');
  if (wind_thick == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('wind_thick = '+wind_thick);
     scao_para=[scao_para wind_thick];
  end

  ground_vel=get(edit323,'String');
  if (ground_vel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('ground_vel = '+ground_vel);
     scao_para=[scao_para ground_vel];
  end

  trop_vel=get(edit324,'String');
  if (trop_vel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('trop_vel = '+trop_vel);
     scao_para=[scao_para trop_vel];
  end

  wf_pixel=get(edit52,'String');
  if (wf_pixel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('wf_pixel = '+wf_pixel);
     scao_para=[scao_para wf_pixel];
  end

  lay_pixel=get(edit53,'String');
  if (lay_pixel == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('lay_pixel = '+lay_pixel);
     scao_para=[scao_para lay_pixel];
  end

  object_n_x=get(edit111,'String');
  if (object_n_x == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('object_n_x = '+object_n_x);
     scao_para=[scao_para object_n_x];
  end

  object_n_y=get(edit112,'String');
  if (object_n_y == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('object_n_y = '+object_n_y);
     scao_para=[scao_para object_n_y];
  end

  object_gap=get(edit113,'String');
  if (object_gap == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('object_gap = '+object_gap);
     scao_para=[scao_para object_gap];
  end

  verbose=get(edit041,'String');
  if (verbose == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('verbose = '+verbose);
     scao_para=[scao_para verbose];
  end

  vverbose=get(edit042,'String');
  if (vverbose == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('vverbose = '+vverbose);
     scao_para=[scao_para vverbose];
  end

  aperture_code=get(edit121,'String');
  if (aperture_code == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('aperture_code = '+aperture_code);
     scao_para=[scao_para aperture_code];
  end

  inner_diameter=get(edit122,'String');
  if (inner_diameter == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('inner_diameter = '+inner_diameter);
     scao_para=[scao_para inner_diameter];
  end

  outer_diameter=get(edit123,'String');
  if (outer_diameter == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('outer_diameter = '+outer_diameter);
     scao_para=[scao_para outer_diameter];
  end

  nspiders=get(edit124,'String');
  if (nspiders == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('nspiders = '+nspiders);
     scao_para=[scao_para nspiders];
  end

  spider_width=get(edit125,'String');
  if (spider_width == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('spider_width = '+spider_width);
     scao_para=[scao_para spider_width];
  end

  edge_length=get(edit126,'String');
  if (edge_length == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('edge_length = '+edge_length);
     scao_para=[scao_para edge_length];
  end

  gap_size=get(edit127,'String');
  if (gap_size == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('gap_size = '+gap_size);
     scao_para=[scao_para gap_size];
  end

  lenslet_focus=get(edit231,'String');
  if (lenslet_focus == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('lenslet_focus = '+lenslet_focus);
     scao_para=[scao_para lenslet_focus];
  end

  lenslet_pitch=get(edit232,'String');
  if (lenslet_pitch == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('lenslet_pitch = '+lenslet_pitch);
     scao_para=[scao_para lenslet_pitch];
  end

  pixels_per_lenslet=get(edit233,'String');
  if (pixels_per_lenslet == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('pixels_per_lenslet = '+pixels_per_lenslet);
     scao_para=[scao_para pixels_per_lenslet];
  end

  pixels_per_xform=get(edit234,'String');
  if (pixels_per_xform == '') then
      x_message(['the parameter aperture_code is not valid']); 
      fin=%t;
  else
     execstr('pixels_per_xform = '+pixels_per_xform);
     scao_para=[scao_para pixels_per_xform];
  end

  geomtry=get(edit031,'String');
  if (geomtry == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('geomtry = '+geomtry);
     scao_para=[scao_para geomtry];
  end

  diffrctive=get(edit032,'String');
  if (diffrctive == '') then
      x_message(['the propagating parameter diffrctive is not valid']); 
      fin=%t;
  else
     execstr('diffrctive = '+diffrctive);
     if ( geomtry + diffrctive <> 1.) then
		x_message(['the propagating method parameter is not accordant']); 
	    fin=%t;
     else
		scao_para=[scao_para diffrctive];
     end
  end
  
  nsubaps=get(edit41,'String');
  if (nsubaps == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('nsubaps = '+nsubaps);
     scao_para=[scao_para nsubaps];
  end

  real_threshold=get(edit42,'String');
  if (real_threshold == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('real_threshold = '+real_threshold);
     scao_para=[scao_para real_threshold];
  end

  eigen_threshold=get(edit43,'String');
  if (eigen_threshold == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('eigen_threshold = '+eigen_threshold);
     scao_para=[scao_para eigen_threshold];
  end

  atm_type=get(edit330,'String');
  if (atm_type == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('atm_type = '+atm_type);
     scao_para=[scao_para atm_type];
  end

  nlayers=get(edit331,'String');
  if (nlayers == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('nlayers = '+nlayers);
     scao_para=[scao_para nlayers];
  end

  r_0=get(edit332,'String');
  if (r_0 == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('r_0 = '+r_0);
     scao_para=[scao_para r_0];
  end

  r_0_wavelength=get(edit333,'String');
  if (r_0_wavelength == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('r_0_wavelength = '+r_0_wavelength);
     scao_para=[scao_para r_0_wavelength];
  end

  pseudowind=get(edit334,'String');
  if (pseudowind == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('pseudowind = '+pseudowind);
     scao_para=[scao_para pseudowind];
  end

  outer_scale=get(edit335,'String');
  if (outer_scale == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('outer_scale = '+outer_scale);
     scao_para=[scao_para outer_scale];
  end

  ground_quality=get(edit336,'String');
  if (ground_quality == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('ground_quality = '+ground_quality);
     scao_para=[scao_para ground_quality];
  end

  anisoplan_quality=get(edit337,'String');
  if (anisoplan_quality == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('anisoplan_quality = '+anisoplan_quality);
     scao_para=[scao_para anisoplan_quality];
  end

  extended_profile=get(edit338,'String');
  if (extended_profile == '') then
      x_message(['the propagating parameter geomtry is not valid']); 
      fin=%t;
  else
     execstr('extended_profile = '+extended_profile);
     scao_para=[scao_para extended_profile];
  end

  for i=1:size(list_heights,'c') 
	execstr('height = '+list_heights(i));
	scao_para=[scao_para height];
  end
  
  for i=1:size(list_objects,'c') 
	execstr('object = '+list_objects(i));
	scao_para=[scao_para object];
  end
  
  if ~fin then
    fprintfMat("scao_para.txt",scao_para,"%e");
    x_message(['all the parameters are valid';
		'now you can execute the simulations '])
  end
  fin=resume(fin);
  
endfunction

///////////////////////////////////////////////////////

function do_execute()

  if fin then
	x_message(['you should check the valid of parameters';
		'please execute the check first '])
  else
	x_message(['start simulating ? ';
		'it will take a while !']);
		
	scao_para=fscanfMat('scao_para.txt');	
	scao_simulation(scao_para);
	
	Dir=pwd();
	x_message(['Now the simulation has finished';
		'all the results has been writen to fits files';
		'the current directory is '+Dir;
		'you can go there to analyse them']);
  end  
  
endfunction

///////////////////////////////////////////////////////


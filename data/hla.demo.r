#$Author: sinnwell $
#
#$Date: 2003/09/08 17:32:55 $
#
#$Header: /people/biostat3/sinnwell/Rdir/Make/RCS/hla.demo.q,v 1.4 2003/09/08 17:32:55 sinnwell Exp $
#
#$Id: hla.demo.q,v 1.4 2003/09/08 17:32:55 sinnwell Exp $
#
#$Locker:  $
#
#$Log: hla.demo.q,v $
#Revision 1.4  2003/09/08 17:32:55  sinnwell
#eleven loci version, with added keywords and license
#
#Revision 1.2  2003/03/06 23:17:21  sinnwell
#add license text
#
#Revision 1.1  2002/09/09 19:53:18  sinnwell
#Initial revision
#
# License: 
# 
# Copyright 2003 Mayo Foundation for Medical Education and Research. 
# 
# This program is free software; you can redistribute it and/or modify it under the terms of 
# the GNU General Public License as published by the Free Software Foundation; either 
# version 2 of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful, but WITHOUT ANY 
# WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for 
# more details.
# 
# You should have received a copy of the GNU General Public License along with this 
# program; if not, write to the Free Software Foundation, Inc., 59 Temple Place, Suite 330, 
# Boston, MA 02111-1307 USA
# 
# For other licensing arrangements, please contact Daniel J. Schaid.
# 
# Daniel J. Schaid, Ph.D.
# Division of Biostatistics
# Harwick Building � Room 775
# Mayo Clinic
# 200 First St., SW
# Rochester, MN 55905
# 
# phone: 507-284-0639
# fax:      507-284-9542
# email: schaid@mayo.edu
# 

"hla.demo"<-
structure(.Data = list(resp = c(3.7010000000000001, 1.179, 1.6599999999999999, 
	1.9670000000000001, 0.79600000000000004, 1.153, 1.3080000000000001, 
	2.0859999999999999, 1.155, 1.5509999999999999, 1.2290000000000001, 
	2.3730000000000002, 2.4300000000000002, 2.2160000000000002, 
	2.3650000000000002, 2.927, 2.5579999999999998, 1.097, 
	2.6339999999999999, 1.7, 1.95, 1.883, 2.3919999999999999, 1.891, 
	2.1970000000000001, 1.399, 2.4430000000000001, 2.1760000000000002, 
	0.63900000000000001, 0.18099999999999999, 0.48299999999999998, 
	0.28000000000000003, 2.0899999999999999, 0.50600000000000001, 
	0.096000000000000002, 0.40200000000000002, 0.35399999999999998, 
	0.63800000000000001, 0.36899999999999999, 0.67400000000000004, 
	1.1799999999999999, 0.70399999999999996, 0.77100000000000002, 
	0.13300000000000001, 0.59599999999999997, 0.78400000000000003, 
	0.70999999999999996, 0.017000000000000001, 0.40100000000000002, 
	0.66400000000000003, 0.67600000000000005, 0.69199999999999995, 
	1.8779999999999999, 0.60699999999999998, 0, 0.043999999999999997, 
	0.72499999999999998, 0.13400000000000001, 0, 0.438, 0.432, 
	0.35799999999999998, 0.75800000000000001, 0.622, 0.36099999999999999, 
	0.71699999999999997, 0.73799999999999999, 0.68500000000000005, 
	0.72899999999999998, 0.76000000000000001, 2.9039999999999999, 
	2.9809999999999999, 3.1200000000000001, 0.68200000000000005, 0.432, 
	0.61499999999999999, 0.53300000000000003, 0.67200000000000004, 
	0.33200000000000002, 0.71899999999999997, 0.33300000000000002, 
	0.61199999999999999, 0.621, 0.378, 0.28899999999999998, 
	0.19900000000000001, 0.68200000000000005, 0.72799999999999998, 
	3.7469999999999999, 3.0310000000000001, 3.266, 0.67100000000000004, 
	0.61099999999999999, 3.2320000000000002, 3.23, 3.7000000000000002, 
	3.601, 3.2709999999999999, 2.5569999999999999, 3.7000000000000002, 
	3.3199999999999998, 3.0609999999999999, 3.1190000000000002, 
	3.0539999999999998, 3.0950000000000002, 2.9289999999999998, 3.423, 
	3.3050000000000002, 3.5750000000000002, 3.379, 3.5670000000000002, 
	2.923, 3.2120000000000002, 3.0830000000000002, 3.0070000000000001, 
	3.431, 3.4449999999999998, 3.1930000000000001, 4.6660000000000004, 
	2.9670000000000001, 3.3340000000000001, 3.226, 3.1869999999999998, 
	2.0699999999999998, 3.274, 3.206, 3.73, 3.1960000000000002, 
	3.1440000000000001, 3.0489999999999999, 2.952, 3.3079999999999998, 
	3.7959999999999998, 2.6890000000000001, 1.605, 3.1949999999999998, 
	1.6020000000000001, 2.9780000000000002, 2.7629999999999999, 0.123, 
	2.1389999999999998, 1.036, 1.117, 2.2160000000000002, 1.958, 1.405, 
	1.1200000000000001, 1.147, 2.5550000000000002, 1.0760000000000001, 
	1.8500000000000001, 1.599, 1.1399999999999999, 2.133, 1.349, 
	1.7729999999999999, 1.9590000000000001, 0.77900000000000003, 
	2.2080000000000002, 4.3289999999999997, 1.0669999999999999, 1.919, 
	1.931, 2.0699999999999998, 2.5099999999999998, 1.387, 
	1.5289999999999999, 1.1739999999999999, 1.923, 1.954, 1.47, 1.893, 
	1.6859999999999999, 2.726, 2.1400000000000001, 0.28999999999999998, 
	1.347, 1.7949999999999999, 1.329, 2.0590000000000002, 1.359, 1.151, 
	1.774, 1.01, 1.044, 1.8140000000000001, 2.097, 1.1020000000000001, 
	3.468, 2.911, 2.9870000000000001, 1.919, 1.208, 1.319, 
	3.4119999999999999, 2.3969999999999998, 4.1040000000000001, 
	3.5499999999999998, 1.119, 1.159, 2.3519999999999999, 
	2.2320000000000002, 1.5449999999999999, 1.2250000000000001, 3.984, 
	1.325, 3.121, 2.9260000000000002, 3.423, 3.569, 3.5750000000000002, 
	3.8599999999999999, 3.133, 3.1059999999999999, 0.41899999999999998, 
	3.6000000000000001, 3.2229999999999999, 0.156, 0.46100000000000002, 
	2.9359999999999999), resp.cat = structure(.Data = c(3, 2, 2, 2, 1, 2, 2,
	2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 
	1, 2, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 1, 1, 1, 1, 
	1, 1, 1, 1, 1, 1, 1, 1, 1, 3, 3, 3, 1, 1, 3, 3, 3, 3, 3, 2, 3, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3, 3, 3, 
	3, 3, 3, 3, 3, 3, 2, 2, 3, 2, 3, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 1, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 3, 3, 2, 2, 2, 3, 2, 3, 3, 2, 
	2, 2, 2, 2, 2, 3, 2, 3, 3, 3, 3, 3, 3, 3, 3, 1, 3, 3, 1, 1, 3), .Label
	 = c("low", "normal", "high"), class = "factor"), male = c(1, 0, 1, 0, 
	0, 1, 0, 1, 1, 1, 0, 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 0, 1, 0, 
	1, 1, 0, 1, 0, 0, 0, 1, 0, 1, 1, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 
	0, 1, 0, 1, 0, 0, 0, 1, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 0, 
	1, 1, 0, 1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 0, 1, 1, 0, 1, 1, 0, 0, 1, 1, 
	0, 0, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1, 1, 0, 1, 1, 1, 1, 0, 1, 1, 
	1, 1, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 1, 0, 0, 1, 0, 
	1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 0, 
	0, 0, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 
	0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 0, 1, 0, 1, 0, 1, 0, 0, 0, 0, 1),
	age = c(15, 15, 17, 21, 15, 15, 15, 15, 15, 20, 15, 14, 15, 15, 15, 15, 
	15, 15, 14, 15, 20, 15, 14, 15, 87, 15, 14, 15, 15, 15, 14, 15, 14, 15, 
	15, 15, 15, 14, 14, 15, 16, 15, 15, 20, 15, 16, 20, 15, 17, 15, 16, 18, 
	15, 18, 17, 15, 16, 15, 15, 17, 22, 16, 15, 16, 15, 15, 15, 15, 15, 24, 
	15, 16, 15, 16, 20, 18, 15, 15, 15, 20, 15, 15, 15, 15, 15, 15, 14, 15, 
	15, 15, 15, 15, 15, 28, 15, 14, 55, 14, 14, 18, 16, 15, 15, 15, 18, 15, 
	15, 15, 15, 19, 15, 14, 16, 15, 16, 15, 16, 15, 15, 15, 14, 15, 15, 15, 
	17, 15, 19, 16, 18, 22, 7, 15, 15, 15, 15, 16, 15, 15, 15, 12, 17, 56, 
	18, 15, 19, 15, 19, 16, 19, 16, 15, 17, 16, 14, 15, 20, 15, 20, 15, 9, 
	18, 14, 14, 19, 15, 16, 15, 22, 15, 16, 17, 15, 15, 15, 16, 15, 17, 18, 
	15, 15, 15, 15, 15, 18, 18, 15, 16, 38, 18, 14, 17, 17, 15, 19, 15, 17, 
	15, 15, 17, 14, 18, 16, 15, 16, 18, 15, 16, 15, 17, 15, 2, 15, 15, 17, 
	15, 15, 16, 15, 14, 15), DPB.a1 = c(2012, 1011, 401, 301, 301, 402, 401,
	401, 2012, 401, 0, 1011, 0, 401, 1011, 1011, 2012, 2012, 401, 2012, 301,
	1011, 1011, 401, 401, 401, 301, 401, 1011, 2012, 401, 1011, 402, 1011, 
	401, 402, 2012, 301, 401, 2102, 401, 401, 401, 301, 402, 2012, 1011, 
	1011, 401, 401, 401, 301, 301, 1011, 2012, 1011, 301, 402, 301, 1011, 
	901, 2012, 1011, 401, 2012, 301, 1011, 401, 2012, 401, 2012, 2012, 402, 
	401, 301, 401, 401, 401, 2012, 1011, 401, 401, 401, 401, 2012, 401, 
	2012, 2012, 1011, 601, 401, 0, 2012, 2012, 301, 301, 301, 1011, 301, 
	301, 401, 301, 2012, 1011, 0, 401, 301, 401, 1011, 1011, 2012, 2012, 
	401, 401, 1011, 2012, 402, 401, 1011, 401, 301, 402, 401, 0, 401, 1011, 
	301, 1011, 1011, 601, 601, 301, 402, 401, 401, 2012, 0, 2012, 401, 1011,
	0, 301, 0, 401, 401, 402, 2012, 0, 401, 301, 401, 1011, 401, 601, 401, 
	301, 1401, 1401, 401, 401, 401, 0, 401, 301, 401, 401, 0, 0, 1011, 301, 
	402, 301, 401, 401, 301, 401, 0, 2012, 401, 0, 0, 402, 401, 202, 401, 
	301, 0, 301, 2012, 2012, 2012, 1011, 0, 0, 401, 2012, 2012, 2012, 5101, 
	402, 301, 202, 2012, 1011, 301, 301, 1011, 401, 301, 401, 401, 0, 501, 
	401, 401, 2012, 1011, 1011, 401, 401), DPB.a2 = c(1301, 401, 401, 11011,
	1001, 11011, 401, 402, 401, 402, 0, 402, 0, 401, 2012, 2301, 401, 401, 
	401, 401, 401, 1011, 402, 401, 601, 402, 1001, 401, 2012, 401, 402, 401,
	2301, 401, 1601, 402, 601, 401, 11011, 401, 401, 1701, 402, 401, 901, 
	2012, 1011, 401, 501, 401, 401, 402, 1401, 301, 1001, 2012, 301, 402, 
	401, 401, 1401, 2012, 401, 402, 401, 11011, 1301, 402, 301, 401, 1901, 
	1301, 601, 1501, 1301, 11011, 401, 1401, 401, 401, 402, 401, 20011, 401,
	401, 401, 401, 1001, 301, 1901, 301, 0, 601, 401, 401, 401, 401, 2012, 
	401, 401, 1601, 401, 401, 2012, 0, 1401, 401, 401, 1701, 301, 401, 401, 
	402, 501, 401, 401, 3301, 401, 301, 501, 401, 1601, 402, 0, 20012, 301, 
	402, 402, 402, 601, 401, 401, 26012, 1501, 401, 301, 0, 301, 401, 301, 
	0, 401, 0, 401, 1001, 402, 20011, 0, 401, 301, 401, 901, 401, 1001, 401,
	401, 1701, 1701, 401, 401, 401, 0, 401, 402, 401, 1401, 0, 0, 401, 401, 
	4801, 11011, 1001, 401, 401, 401, 0, 402, 401, 0, 0, 402, 401, 401, 401,
	4701, 0, 401, 401, 301, 1301, 401, 0, 0, 401, 2012, 402, 402, 5501, 
	1001, 402, 401, 402, 2012, 401, 402, 401, 601, 4001, 402, 402, 0, 1301, 
	1401, 301, 401, 401, 301, 401, 401), DPA.a1 = c(103, 103, 103, 103, 103,
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 201, 103, 103, 103, 103, 201, 103, 103, 103, 103, 103, 0, 103,
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 201, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 201, 103, 
	103, 103, 103, 103, 201, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 201, 103, 103, 201, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 0, 103, 103, 103, 103, 103, 103, 103, 103, 103,
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 201, 103, 103, 103, 103, 103, 201, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 103, 
	103, 103, 103, 103, 103, 103, 103, 103, 103, 202, 103, 103, 103, 103, 
	103, 103, 103), DPA.a2 = c(201, 201, 103, 201, 201, 201, 103, 103, 103, 
	103, 103, 201, 202, 103, 201, 201, 103, 103, 103, 103, 103, 201, 201, 
	103, 103, 103, 202, 103, 201, 103, 103, 201, 0, 201, 103, 103, 103, 103,
	201, 103, 103, 201, 103, 103, 201, 103, 201, 201, 202, 103, 103, 103, 
	201, 201, 201, 201, 103, 103, 103, 201, 202, 103, 201, 103, 103, 201, 
	201, 103, 103, 103, 202, 201, 103, 104, 103, 201, 103, 201, 103, 201, 
	103, 103, 103, 103, 103, 103, 103, 201, 202, 202, 201, 103, 103, 103, 
	103, 103, 103, 201, 103, 103, 103, 103, 103, 201, 103, 202, 103, 103, 
	201, 201, 103, 103, 103, 202, 201, 103, 103, 103, 201, 202, 103, 103, 
	103, 0, 103, 201, 103, 201, 201, 103, 103, 103, 103, 104, 103, 103, 103,
	103, 103, 201, 201, 103, 201, 103, 103, 103, 103, 103, 103, 103, 103, 
	201, 103, 201, 103, 103, 103, 201, 103, 103, 103, 103, 103, 103, 103, 
	201, 103, 103, 201, 103, 103, 201, 201, 103, 103, 103, 103, 201, 103, 
	103, 201, 103, 103, 103, 103, 103, 103, 103, 103, 103, 201, 201, 103, 
	103, 103, 103, 103, 103, 201, 201, 103, 103, 103, 201, 103, 103, 202, 
	103, 301, 103, 103, 103, 202, 201, 103, 103, 201, 201, 103, 103), 
	DMA.a1 = c(101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 0, 101, 101, 101, 101, 101, 101, 101, 101,
	102, 101, 101, 101, 104, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 102, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 104, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101), DMA.a2 = c(101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 102, 103, 102, 101, 101, 101, 103, 101, 101, 101, 101, 101, 
	104, 101, 101, 0, 101, 101, 101, 102, 101, 101, 101, 101, 104, 101, 101,
	103, 104, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 104, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 102, 101, 
	103, 101, 101, 101, 101, 102, 102, 101, 101, 101, 102, 103, 101, 101, 
	101, 104, 102, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 102, 
	101, 102, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 102, 
	101, 104, 101, 103, 101, 103, 101, 102, 102, 102, 103, 101, 103, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 104, 
	101, 101, 102, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 104, 
	101, 102, 103, 102, 101, 101, 101, 101, 101, 102, 101, 101, 101, 101, 
	101, 101, 101, 102, 101, 101, 101, 101, 101, 101, 102, 101, 102, 101, 
	101, 101, 101, 102, 101, 104, 101, 101, 101, 101, 101, 102, 101, 101, 
	101, 101, 101, 101, 101, 104, 101, 101, 102, 103, 102, 101, 101, 102, 
	101, 103, 101, 101, 101, 101, 101, 101), DMB.a1 = c(101, 101, 101, 101, 
	101, 101, 101, 101, 103, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 0,
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 103, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 103, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 0, 101, 101, 101, 101, 101, 101, 101, 103,
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 103, 103, 101, 101, 101, 
	101, 101, 101, 101, 103, 101, 101, 101, 101, 101, 101, 103, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 0, 101, 101, 101, 101, 101,
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101), DMB.a2 = c(101, 101, 101, 104, 101, 101, 101, 101, 103, 
	101, 103, 101, 101, 103, 101, 101, 101, 101, 104, 101, 101, 101, 101, 
	101, 101, 101, 101, 101, 101, 101, 101, 101, 0, 101, 101, 101, 101, 101,
	102, 101, 101, 101, 101, 101, 104, 101, 101, 104, 104, 103, 103, 101, 
	101, 101, 102, 101, 101, 101, 101, 102, 101, 101, 101, 103, 101, 104, 
	101, 101, 101, 103, 101, 102, 101, 101, 103, 101, 101, 101, 101, 101, 
	101, 103, 101, 103, 103, 102, 101, 101, 101, 101, 101, 101, 101, 101, 
	101, 101, 101, 101, 103, 101, 103, 103, 101, 101, 101, 101, 104, 103, 
	101, 101, 101, 101, 101, 103, 101, 101, 102, 101, 101, 101, 103, 101, 
	101, 0, 101, 101, 103, 101, 101, 101, 101, 103, 101, 101, 103, 102, 101,
	103, 101, 101, 101, 101, 101, 103, 101, 101, 101, 101, 101, 101, 103, 
	101, 101, 101, 101, 103, 103, 101, 101, 102, 101, 101, 101, 101, 103, 
	101, 101, 101, 103, 101, 101, 104, 101, 101, 101, 103, 101, 101, 101, 
	101, 101, 101, 103, 0, 103, 103, 101, 101, 103, 103, 103, 101, 101, 103,
	101, 101, 101, 101, 102, 101, 101, 101, 101, 101, 101, 101, 102, 101, 
	101, 101, 102, 101, 101, 101, 101, 103, 101, 101, 101, 101), TAP1.a1 = 
	structure(.Data = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 
	3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2), .Label = c("0", "A", "B"), class = 
	"factor"), TAP1.a2 = structure(.Data = c(2, 2, 2, 3, 3, 2, 2, 2, 2, 2, 
	2, 4, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 3, 3, 3, 3, 2, 2, 2, 2, 2, 1, 3, 
	2, 3, 2, 2, 2, 2, 2, 2, 4, 2, 3, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 3, 3, 
	2, 2, 2, 2, 2, 3, 3, 3, 2, 4, 3, 2, 4, 2, 3, 4, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 3, 2, 3, 2, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 2, 
	3, 2, 2, 2, 2, 4, 3, 2, 2, 2, 2, 3, 2, 3, 3, 2, 2, 2, 4, 3, 2, 2, 2, 4, 
	2, 2, 2, 3, 2, 3, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 3, 2, 2, 2, 3, 3, 2, 4, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 
	3, 2, 2, 3, 2, 2, 4, 3, 2, 3, 2, 2, 3, 2, 2, 2, 3, 3, 4, 2, 2, 2, 4, 2, 
	2, 1, 2, 2, 3, 2, 2, 2, 2, 3, 2, 3, 2, 2, 3, 3, 3, 2), .Label = c("0", 
	"A", "B", "C"), class = "factor"), TAP2.a1 = structure(.Data = c(2, 2, 
	2, 2, 2, 2, 1, 2, 2, 3, 1, 2, 2, 2, 2, 2, 2, 2, 1, 2, 1, 2, 2, 1, 2, 2, 
	2, 2, 2, 2, 1, 2, 1, 2, 2, 2, 2, 3, 2, 2, 4, 3, 2, 2, 1, 2, 2, 2, 2, 4, 
	3, 2, 3, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 4, 2, 2, 
	2, 2, 2, 2, 2, 2, 3, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 
	2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 4, 3, 2, 2, 3, 2, 2, 1, 2, 2, 2, 2, 
	1, 1, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 3, 2, 2, 2, 2, 2, 3, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 1, 2, 2, 2, 1, 2, 2, 3, 1, 2, 1, 2, 
	2, 2, 2, 2, 2, 2, 2, 2, 1, 3, 1, 2, 1, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 
	2, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 1, 2, 2, 1, 2, 
	2, 2), .Label = c("0", "A", "B", "C"), class = "factor"), TAP2.a2 = 
	structure(.Data = c(3, 2, 3, 2, 2, 2, 1, 3, 3, 3, 1, 2, 3, 2, 3, 3, 3, 
	6, 1, 4, 1, 2, 3, 1, 2, 3, 3, 2, 3, 2, 1, 3, 1, 2, 3, 3, 3, 3, 5, 3, 4, 
	3, 2, 2, 1, 3, 2, 2, 3, 5, 3, 2, 3, 3, 4, 2, 3, 2, 2, 2, 2, 2, 6, 4, 3, 
	2, 3, 2, 4, 1, 2, 4, 4, 2, 4, 2, 2, 4, 2, 2, 3, 3, 3, 1, 3, 2, 2, 3, 4, 
	3, 2, 3, 6, 4, 1, 3, 2, 3, 3, 2, 1, 3, 2, 3, 4, 3, 2, 2, 2, 2, 2, 3, 3, 
	2, 3, 3, 3, 1, 2, 4, 2, 3, 1, 1, 2, 2, 3, 2, 3, 1, 4, 2, 2, 2, 2, 2, 4, 
	2, 3, 3, 3, 2, 2, 3, 3, 3, 2, 4, 2, 3, 3, 3, 3, 2, 3, 2, 1, 3, 1, 3, 2, 
	2, 1, 2, 2, 3, 1, 3, 1, 3, 2, 2, 2, 4, 4, 3, 3, 6, 1, 3, 1, 2, 1, 3, 1, 
	3, 2, 3, 3, 3, 3, 2, 4, 4, 4, 2, 2, 3, 3, 2, 3, 1, 2, 3, 3, 3, 3, 2, 3, 
	2, 3, 3, 3, 1, 2, 2, 1, 3, 4, 3), .Label = c("0", "A", "B", "C", "D", 
	"E"), class = "factor"), DQB.a1 = c(31, 21, 31, 21, 31, 21, 62, 61, 62, 
	51, 31, 62, 63, 62, 51, 51, 31, 62, 62, 31, 31, 21, 51, 42, 64, 51, 63, 
	64, 31, 31, 31, 63, 51, 21, 31, 63, 51, 63, 31, 63, 62, 21, 21, 32, 31, 
	21, 21, 51, 62, 53, 51, 64, 31, 51, 31, 63, 63, 63, 51, 21, 61, 21, 32, 
	21, 31, 63, 21, 21, 53, 21, 21, 53, 62, 31, 52, 21, 51, 31, 21, 31, 51, 
	21, 51, 31, 21, 21, 62, 31, 21, 63, 21, 21, 21, 21, 21, 51, 21, 63, 51, 
	32, 31, 63, 62, 64, 32, 21, 21, 32, 21, 21, 31, 31, 51, 33, 51, 33, 62, 
	31, 32, 51, 21, 21, 51, 31, 32, 33, 21, 62, 21, 51, 21, 62, 21, 31, 31, 
	32, 0, 21, 51, 21, 21, 51, 21, 51, 51, 51, 32, 61, 21, 31, 61, 51, 63, 
	31, 51, 61, 51, 51, 32, 61, 62, 52, 32, 31, 62, 61, 62, 31, 21, 51, 21, 
	21, 51, 31, 62, 21, 51, 21, 51, 51, 21, 31, 31, 21, 31, 31, 21, 63, 31, 
	31, 21, 21, 21, 31, 31, 21, 21, 63, 51, 31, 63, 31, 21, 63, 32, 51, 51, 
	21, 31, 63, 21, 64, 51, 63, 21, 63, 31, 51, 53, 62), DQB.a2 = c(32, 62, 
	63, 31, 42, 21, 64, 31, 63, 51, 32, 32, 63, 32, 63, 62, 63, 31, 42, 64, 
	52, 21, 21, 42, 32, 64, 31, 62, 21, 62, 33, 32, 62, 62, 31, 32, 62, 31, 
	31, 21, 62, 33, 21, 62, 32, 62, 21, 21, 21, 31, 21, 21, 31, 21, 31, 21, 
	62, 32, 21, 32, 61, 31, 62, 31, 64, 21, 33, 32, 33, 31, 64, 62, 53, 62, 
	31, 21, 64, 31, 21, 21, 63, 33, 51, 31, 63, 32, 31, 31, 31, 32, 32, 21, 
	33, 31, 32, 32, 32, 21, 31, 63, 62, 33, 64, 62, 62, 21, 21, 62, 53, 62, 
	51, 62, 51, 62, 33, 62, 51, 62, 62, 51, 64, 32, 62, 32, 32, 64, 31, 21, 
	32, 31, 31, 63, 53, 31, 62, 32, 0, 64, 53, 31, 64, 63, 31, 21, 33, 32, 
	53, 31, 62, 32, 31, 21, 31, 32, 31, 31, 31, 32, 32, 61, 32, 52, 32, 31, 
	62, 31, 62, 32, 31, 21, 21, 21, 51, 32, 62, 51, 63, 31, 52, 63, 32, 64, 
	62, 21, 32, 64, 21, 63, 32, 32, 64, 32, 32, 31, 52, 32, 32, 31, 32, 32, 
	32, 62, 52, 21, 51, 62, 51, 32, 62, 63, 21, 62, 63, 52, 63, 62, 64, 21, 
	21, 51), DQA.a1 = c(301, 102, 103, 201, 401, 201, 102, 501, 102, 101, 
	301, 102, 103, 102, 101, 101, 103, 102, 102, 0, 102, 501, 101, 301, 102,
	101, 103, 102, 501, 102, 201, 103, 101, 401, 501, 103, 102, 103, 301, 
	103, 102, 201, 501, 102, 301, 102, 501, 101, 102, 104, 101, 102, 501, 
	101, 501, 103, 102, 102, 101, 301, 103, 501, 102, 501, 102, 102, 201, 
	301, 104, 201, 102, 102, 102, 102, 102, 201, 101, 301, 201, 301, 101, 
	301, 101, 301, 102, 201, 501, 501, 201, 103, 301, 201, 201, 501, 301, 
	104, 102, 501, 101, 103, 102, 103, 102, 102, 102, 401, 201, 102, 104, 
	102, 101, 103, 101, 102, 104, 102, 101, 104, 102, 104, 102, 201, 102, 0,
	102, 102, 201, 102, 101, 101, 501, 102, 104, 401, 102, 301, 104, 102, 
	104, 501, 102, 102, 401, 0, 101, 101, 104, 102, 0, 301, 103, 101, 103, 
	401, 101, 102, 101, 101, 101, 103, 102, 102, 301, 102, 102, 102, 102, 
	201, 102, 101, 201, 201, 101, 301, 102, 104, 101, 201, 101, 101, 201, 
	102, 102, 501, 301, 102, 301, 103, 301, 201, 102, 301, 301, 501, 401, 
	301, 201, 103, 101, 301, 103, 102, 102, 103, 101, 101, 101, 201, 102, 
	102, 102, 102, 103, 102, 102, 102, 102, 101, 104, 101), DQA.a2 = c(501, 
	201, 501, 501, 501, 501, 501, 501, 103, 101, 301, 301, 501, 301, 103, 
	103, 501, 501, 401, 0, 301, 501, 501, 401, 301, 102, 501, 102, 501, 501,
	501, 301, 102, 501, 501, 301, 104, 301, 301, 501, 102, 201, 501, 301, 
	501, 201, 501, 501, 501, 301, 201, 501, 601, 201, 501, 501, 401, 301, 
	201, 501, 301, 501, 301, 501, 501, 201, 501, 501, 301, 301, 201, 104, 
	104, 301, 601, 501, 102, 301, 501, 501, 103, 501, 101, 501, 501, 301, 
	501, 501, 501, 301, 501, 501, 301, 501, 501, 301, 301, 102, 501, 301, 
	301, 301, 102, 102, 301, 501, 501, 301, 201, 501, 501, 501, 401, 301, 
	201, 301, 102, 301, 301, 104, 501, 301, 104, 0, 301, 301, 501, 501, 501,
	501, 501, 103, 501, 501, 501, 501, 501, 501, 101, 501, 501, 103, 501, 0,
	201, 301, 301, 501, 0, 601, 501, 501, 501, 501, 501, 501, 301, 201, 301,
	103, 301, 301, 401, 501, 201, 301, 102, 501, 501, 501, 501, 501, 101, 
	501, 102, 501, 103, 501, 103, 103, 501, 501, 501, 501, 301, 501, 501, 
	401, 501, 501, 201, 501, 501, 501, 501, 501, 301, 501, 201, 301, 301, 
	301, 501, 501, 301, 102, 104, 301, 301, 102, 501, 102, 104, 103, 501, 
	103, 301, 501, 501, 102), DRB.a1 = c(4, 2, 1, 7, 8, 3, 2, 13, 2, 1, 4, 
	2, 13, 2, 1, 1, 7, 2, 2, 4, 2, 3, 1, 4, 2, 1, 11, 2, 3, 2, 7, 4, 1, 3, 
	11, 4, 2, 4, 4, 3, 2, 7, 3, 2, 4, 2, 3, 1, 2, 4, 1, 3, 8, 1, 11, 3, 2, 
	2, 1, 4, 2, 3, 2, 3, 13, 7, 7, 3, 9, 4, 7, 2, 2, 2, 2, 3, 1, 4, 7, 4, 1,
	9, 1, 4, 2, 4, 11, 11, 7, 4, 3, 3, 7, 3, 3, 4, 2, 2, 1, 4, 2, 13, 2, 11,
	2, 3, 3, 2, 7, 2, 1, 11, 1, 2, 7, 2, 1, 4, 2, 1, 3, 4, 2, 4, 2, 10, 7, 
	2, 1, 1, 3, 2, 8, 8, 2, 4, 11, 7, 1, 3, 3, 2, 3, 3, 1, 1, 4, 2, 2, 4, 
	11, 1, 1, 11, 1, 8, 1, 1, 1, 2, 2, 2, 8, 2, 2, 2, 2, 7, 8, 1, 3, 7, 1, 
	8, 2, 3, 1, 4, 1, 1, 3, 11, 8, 8, 4, 13, 4, 8, 4, 7, 7, 4, 3, 3, 8, 4, 
	4, 11, 1, 4, 4, 2, 2, 10, 1, 1, 1, 4, 4, 2, 2, 2, 10, 2, 2, 2, 4, 1, 3, 
	1), DRB.a2 = c(11, 7, 13, 7, 11, 7, 13, 11, 2, 1, 4, 4, 13, 4, 13, 13, 
	13, 14, 8, 13, 4, 3, 3, 8, 4, 13, 13, 13, 11, 11, 14, 13, 2, 8, 11, 13, 
	10, 13, 4, 13, 2, 7, 3, 4, 8, 7, 3, 3, 3, 14, 7, 13, 13, 7, 11, 13, 8, 
	4, 7, 3, 4, 11, 4, 11, 13, 13, 3, 4, 14, 7, 13, 9, 14, 4, 2, 7, 13, 4, 
	3, 3, 13, 9, 1, 3, 3, 7, 2, 11, 11, 13, 4, 7, 9, 11, 4, 10, 4, 3, 11, 
	13, 4, 9, 13, 13, 9, 7, 7, 4, 14, 3, 11, 4, 8, 9, 10, 9, 2, 14, 4, 14, 
	13, 7, 10, 4, 4, 13, 11, 3, 3, 13, 11, 13, 3, 8, 11, 3, 14, 13, 14, 13, 
	13, 2, 11, 4, 7, 4, 14, 11, 3, 8, 13, 1, 13, 4, 11, 11, 4, 7, 4, 2, 4, 
	4, 4, 4, 7, 4, 10, 8, 13, 3, 7, 7, 1, 11, 2, 8, 4, 7, 2, 13, 7, 11, 8, 
	3, 10, 13, 4, 13, 11, 11, 13, 3, 3, 11, 11, 3, 7, 13, 7, 4, 13, 4, 3, 
	13, 4, 2, 10, 7, 11, 10, 3, 13, 13, 13, 3, 13, 13, 3, 14, 2), B.a1 = c(
	62, 7, 27, 7, 51, 8, 7, 44, 7, 27, 44, 55, 44, 7, 51, 51, 55, 7, 7, 60, 
	44, 8, 8, 35, 51, 7, 62, 57, 8, 7, 7, 38, 18, 8, 35, 60, 7, 44, 44, 44, 
	7, 13, 8, 44, 62, 44, 8, 8, 8, 44, 44, 8, 35, 13, 18, 8, 18, 44, 62, 8, 
	52, 44, 60, 8, 7, 44, 8, 8, 35, 44, 7, 51, 27, 35, 44, 8, 8, 44, 8, 8, 
	0, 51, 27, 8, 7, 7, 44, 49, 44, 62, 7, 8, 50, 8, 7, 7, 18, 7, 51, 8, 7, 
	7, 7, 55, 51, 7, 7, 44, 44, 7, 51, 7, 7, 7, 57, 44, 7, 51, 7, 44, 8, 62,
	18, 13, 44, 51, 57, 8, 8, 35, 7, 7, 8, 45, 44, 7, 0, 44, 38, 8, 8, 51, 
	51, 27, 7, 35, 51, 27, 7, 52, 38, 8, 44, 44, 44, 44, 44, 39, 62, 62, 7, 
	39, 7, 14, 7, 44, 51, 7, 8, 8, 45, 7, 8, 7, 8, 7, 44, 7, 62, 7, 57, 51, 
	39, 8, 27, 8, 62, 44, 62, 57, 62, 8, 51, 8, 35, 14, 14, 51, 44, 7, 7, 
	13, 62, 8, 18, 51, 27, 44, 45, 7, 7, 7, 7, 38, 7, 7, 7, 8, 8, 7), B.a2
	 = c(61, 44, 62, 44, 55, 44, 63, 44, 35, 35, 60, 27, 62, 60, 7, 35, 61, 
	63, 18, 61, 27, 44, 60, 60, 7, 35, 37, 35, 27, 49, 58, 60, 35, 35, 37, 
	55, 37, 44, 44, 70, 44, 57, 8, 45, 60, 18, 8, 58, 18, 55, 13, 44, 41, 
	44, 56, 60, 18, 62, 47, 60, 61, 35, 60, 60, 14, 60, 57, 8, 60, 13, 60, 
	18, 55, 35, 56, 44, 14, 35, 35, 44, 0, 8, 35, 44, 8, 44, 44, 38, 62, 60,
	8, 13, 60, 35, 8, 62, 60, 8, 62, 44, 44, 48, 63, 35, 7, 8, 18, 62, 38, 
	8, 27, 35, 39, 46, 37, 60, 35, 60, 62, 56, 51, 50, 18, 51, 7, 7, 38, 8, 
	35, 41, 8, 7, 57, 60, 62, 35, 0, 35, 35, 57, 35, 7, 8, 35, 35, 60, 60, 
	60, 8, 60, 39, 35, 27, 60, 27, 35, 27, 14, 7, 61, 35, 56, 62, 57, 57, 
	44, 35, 44, 63, 35, 18, 62, 44, 35, 18, 8, 48, 41, 27, 44, 57, 38, 37, 
	18, 39, 14, 62, 58, 61, 18, 63, 44, 8, 35, 60, 49, 18, 7, 39, 60, 35, 
	60, 35, 44, 60, 7, 14, 60, 37, 44, 8, 60, 35, 42, 8, 38, 27, 48, 44, 35
	), A.a1 = c(11, 3, 2, 3, 2, 1, 1, 24, 3, 2, 2, 2, 3, 2, 3, 2, 2, 3, 3, 
	1, 2, 1, 2, 2, 3, 1, 1, 2, 1, 2, 2, 2, 3, 1, 1, 28, 2, 2, 2, 2, 2, 2, 1,
	3, 2, 3, 1, 1, 1, 2, 1, 1, 2, 2, 32, 1, 2, 3, 2, 1, 11, 2, 2, 1, 28, 2, 
	1, 1, 2, 2, 2, 2, 2, 3, 2, 1, 3, 2, 1, 2, 0, 1, 2, 1, 1, 1, 1, 26, 29, 
	1, 2, 1, 2, 1, 1, 2, 2, 1, 2, 1, 11, 3, 32, 23, 2, 2, 2, 2, 2, 1, 3, 1, 
	2, 11, 2, 2, 29, 2, 2, 2, 1, 3, 2, 30, 32, 11, 1, 1, 1, 3, 1, 3, 1, 3, 
	32, 3, 0, 2, 2, 2, 1, 2, 2, 2, 3, 3, 3, 11, 1, 2, 2, 1, 2, 29, 2, 2, 2, 
	32, 2, 3, 2, 24, 2, 1, 1, 2, 11, 2, 11, 1, 2, 1, 1, 11, 1, 1, 1, 2, 2, 
	2, 1, 3, 2, 1, 2, 1, 3, 2, 1, 2, 11, 1, 1, 1, 2, 30, 25, 3, 2, 2, 3, 1, 
	2, 1, 2, 3, 31, 31, 1, 3, 2, 2, 29, 28, 1, 3, 2, 1, 2, 3), A.a2 = c(24, 
	23, 2, 30, 2, 2, 3, 30, 3, 11, 2, 23, 24, 3, 11, 3, 26, 24, 24, 2, 32, 
	3, 26, 31, 28, 3, 24, 26, 31, 26, 3, 28, 25, 3, 2, 31, 3, 28, 2, 23, 25,
	26, 2, 24, 3, 25, 1, 1, 32, 11, 32, 2, 26, 30, 24, 2, 25, 29, 3, 1, 24, 
	24, 32, 2, 24, 2, 2, 11, 3, 24, 24, 25, 2, 11, 2, 29, 26, 32, 11, 2, 0, 
	1, 3, 2, 3, 2, 26, 26, 24, 2, 24, 30, 31, 28, 29, 2, 25, 1, 31, 1, 24, 
	24, 33, 24, 3, 2, 30, 3, 31, 2, 31, 24, 3, 26, 26, 32, 31, 28, 3, 24, 1,
	30, 3, 32, 24, 31, 2, 1, 3, 26, 2, 3, 32, 29, 24, 24, 0, 3, 26, 29, 28, 
	24, 2, 3, 11, 31, 11, 11, 32, 2, 26, 28, 11, 31, 3, 11, 28, 33, 2, 28, 
	11, 24, 3, 3, 11, 31, 24, 2, 32, 11, 30, 31, 2, 32, 25, 23, 26, 23, 3, 
	3, 2, 26, 24, 2, 32, 28, 29, 29, 2, 24, 24, 33, 28, 32, 31, 33, 32, 29, 
	3, 2, 11, 2, 24, 2, 3, 3, 23, 23, 2, 11, 30, 32, 24, 24, 3, 3, 31, 32, 
	24, 3)), row.names = c("1", "2", "3", "4", "5", "6", "7", "9", "11", 
	"12", "14", "15", "16", "17", "18", "19", "20", "21", "22", "23", "24", 
	"25", "26", "27", "28", "29", "30", "31", "34", "35", "36", "37", "38", 
	"39", "40", "41", "42", "43", "44", "46", "48", "50", "51", "52", "53", 
	"54", "55", "56", "57", "58", "59", "60", "62", "63", "64", "65", "66", 
	"67", "68", "69", "70", "71", "72", "73", "74", "75", "76", "77", "78", 
	"79", "80", "81", "82", "83", "84", "85", "87", "88", "89", "91", "93", 
	"94", "95", "96", "97", "98", "99", "100", "101", "102", "103", "104", 
	"105", "106", "107", "108", "109", "110", "111", "113", "114", "115", 
	"116", "118", "119", "120", "121", "122", "123", "125", "126", "127", 
	"129", "130", "131", "132", "133", "134", "136", "137", "138", "139", 
	"140", "141", "142", "143", "144", "145", "146", "147", "148", "149", 
	"150", "151", "152", "153", "154", "155", "156", "157", "158", "159", 
	"160", "161", "162", "163", "164", "165", "166", "167", "168", "169", 
	"170", "171", "172", "173", "174", "175", "176", "177", "178", "179", 
	"180", "181", "182", "183", "184", "185", "186", "187", "189", "190", 
	"191", "193", "194", "195", "196", "197", "199", "200", "201", "202", 
	"203", "204", "205", "206", "207", "208", "209", "210", "211", "213", 
	"214", "215", "216", "218", "219", "220", "221", "222", "223", "224", 
	"225", "226", "227", "228", "229", "230", "231", "232", "233", "234", 
	"235", "236", "237", "238", "239", "240", "241", "242"), class = 
	"data.frame")

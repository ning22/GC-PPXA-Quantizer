    Copyright 2012 Anna Jezierska
	------------------------------------------------------------------------------------
	GC-PPXA-QUANTIZER is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

	------------------------------------------------------------------------------------
	
	This software implements the vector quantization algorithm
	described in
	
	[1] "A Spatial Regularization Approach for Vector Quantization",
	 C. Chaux, A. Jezierska, J.-C. Pesquet, and H. Talbot,
	 Journal of Mathematical Imaging and Vision, vol. 41, pp. 23-38, 2011
	
	[2] "Image quantization under spatial smoothness constraints" ,
	 A. Jezierska,C. Chaux,  J.-C. Pesquet, and H. Talbot, 
	 International Conference on Image Processing (ICIP), Honk Kong, 26-29 September 2010. 
	
	This software was developed by Anna Jezierska version of 12.2011
	(anna.jezierska@univ-paris-est.fr).
	
	This software uses Graph-Cut library provided by Y. Boykov and V. Kolmogorov, which is described in:
	
	[3]"An Experimental Comparison of Min-Cut/Max-Flow Algorithms for Energy Minimization in Vision."
	Y. Boykov and V. Kolmogorov, In IEEE Transactions on Pattern Analysis and Machine Intelligence (PAMI),
	September 2004
	
	and available at V. Kolmogorov website: http://pub.ist.ac.at/~vnk/software.html (toolbox called "maxflow-v3.01.src.tar.gz").
	The downloaded files: block.h, graph.h, graph.cpp, maxflow.cpp, instances.inc should be saved 
	in the same folder as files of the toolbox provided here.
	
	The color image quantization is initialized with median-cut method introduced in:
	
	[4] P. Heckbert. Color image quantization for frame buffer display.
    SIGGRAPH Comput. Graph., 16(3):297{307, Jul. 1982.
	
	If you use this software for research purposes, you should cite
	the aforementioned paper in any resulting publication.
	
	------------------------------------------------------------------------------------
	SOFTWARE DESCRIPTION:
	This software allows to quantize greyscale (1 channel) and color (3 channels) images.
	The considered criterion is composed with 
	-> data fidelity term
	% (1) mean absolute value criterion
	% (2) quadratic norm
	-> a regularization term favorizing piecewise constant images
	rho(x) = weight sum(phi(psi(x)))
	where psi is a linear operator computing the first-order differences between pixels
	(i.e. the discrete gradient computed in the horizontal and
	vertical directions) and phi is one of the following potential functions
	
	%(1) phi(u) = sqrt(u^2);                %CONVEX 
	%(2) phi(u) = u^2;                      %CONVEX
	%(3) Potts model : 
		 phi(u) = 
	               1 if  u = 0,  
				   0 otherwise              %NON CONVEX 
 
	The quantizations process uses the PPXA+ algorithm and
    %(1) Graph-Cut alpha expansion in non-convex case
    %(2) Graph-Cut gradient descent algorithm in convex case
    ------------------------------------------------------------------------------------
	USAGE DESCRIPTION:
	
	Requirements:
	 - The GNU Scientific Library (GSL) available at http://www.gnu.org/s/gsl/.
	 - g++ compiler
	 
	Provide following changes into attached Makefile:
     - set path to libraries libgsl.a  and libgslcblas.a by settings GSLLIB and GSLCBLASLIB
	 
	Launching a program:
	1. Compile using make
	2. Start program with following arguments: 
	 <size x> <size y> <nChannels> <ch1> <ch2> <ch3> <Q> <data fidelity> <reg term> <weight> <niter>
	     size x, size y - width and high of image, positive integer 
		 nChannels - number of channels, either 1 or 3
		 ch1, ch2, ch3 -image data, .txt format
	     Q - number of codevectors, positive integer 
	     data fidelity - 1 -> L1, 2-> (L2)^2 
	     reg term (anisotropic TV) - 1 -> L1, 2-> (L2)^2 , 3 -> Potts
	     weight - positive double 
	     niter - max number of iterations ,positive integer
	 
	Example:
    (1)./program 512 512 1 image.txt NULL NULL 40 1 1 20 10 
	(2)./program 512 512 3 ch1.txt ch2.txt ch3.txt 40 1 1 20 10 
	
	Examples (1,2) illustrate quantization of 1-channel and quantization of 3-channel image of size 512x512, respectively.
	
	The resulting images are saved as:
	
	- out_ch1.txt in case of 1-channel image quantization
	- out_ch1.txt, out_ch2.txt ,out_ch3.txt in case of 3-channel image quantization
	------------------------------------------------------------------------------------

	
	

	
	
	
	
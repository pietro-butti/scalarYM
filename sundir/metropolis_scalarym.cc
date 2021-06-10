	//true--->accept, false---reject
bool boltzmann_coin_flip(double deltaS) { 
	if (deltaS<0) return true;
	else{
		double coin = xx.rand();

		if (coin<exp(-deltaS)) return true;
		else return false;
	}
}


double Metropolis_for_one_link(int dir, int site, double tuner, dc* ufield_new) {
	// This is an ULTRA LOCAL Metropolis change on the gauge variables, given one link
	// Generate a random SU(2) matrix X=1+tuner*R,
	// Generate proposal as U(mu,x) ---> X.U(mu,x)
	// 

	
	// generate SU(2) random matrix X=1-tuner*R,
	dc X[4], X_minus_one[4];
	X[0] = 1. + tuner*dc(xx.rand(.5),xx.rand(.5)-1.);
	X[1] =      tuner*dc(xx.rand(.5),xx.rand(.5)-1.);
	X[2] = dc(0.,0.);
	X[3] = dc(0.,0.); 
	norm_su2(X);

    X_minus_one[0] = X[0] - 1.;
    X_minus_one[1] = X[1]     ;
    X_minus_one[2] = X[2]     ;
    X_minus_one[3] = X[3] - 1.;


	// Generate proposal
	dc Uold[4], Unew[4];
	for(int i=0; i<4; i++) Uold[i] = ufield[(dir*nsites+site)*Ncolsquare+i];
	mult_C_equals_AB_for_SU2(Unew,X,Uold);


	// Calculate

	/*
	// Calculate deltaS = beta/2N*tr[ (X-1)*(U(mu,n) + sum_(other i) U(mu,n+i)) ]
	dc UU[4], SS[4];
	for(int i=0; i<4; i++) {
		UU[i] = ufield[(dir*nsites+site)*Ncolsquare+i];

		for(int mu=0; mu<dim; mu++){
			if(mu!=dir) {
				int n_minus_i = neighbor_minus[mu*nsites+site];
				UU[i] += ufield[(dir*nsites+n_minus_i)*Ncolsquare+i];
			} 
		}
	} 
	mult_C_equals_AB_for_SU2(SS,X_minus_one,UU);
	double deltaS = real((SS[0] + SS[3]))*beta/4.;
	*/

	// Update exiting proposal
	for(int i=0; i<4; i++) ufield_new[(dir*nsites+site)*Ncolsquare+i] = Unew[i];

	return deltaS;
}


// This performs a sweep (== update every link variable once) and then perform metropolis
bool Metropolis(double tuner) {

	// Create a bin and fill with gauge variables
	dc *proposal = new dc[ufielddimension];
	for(int ii=0; ii<ufielddimension; ii++) proposal[ii] = ufield[ii];
	
	// Perform the sweep
	double deltaS = 0.;
	for(int site=0; site<nsites; site++)
	for(int dir=0; dir<dim; dir++) {
		deltaS += Metropolis_for_one_link(site,dir,tuner,proposal);
	}


	// Perform Metropolis
	bool coin = boltzmann_coin_flip(deltaS);

	if (coin==true) { //accept
		for(int ii=0; ii<ufielddimension; ii++)  ufield[ii] = proposal[ii];
		return true;
	}
	else return false; //reject

}

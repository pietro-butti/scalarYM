	//true--->accept, false---reject
bool boltzmann_coin_flip(double deltaS) { 
	double coin = xx.rand();
	double probability = min(1.,exp(-deltaS));

	if (coin<probability) return true;
	else return false;
	/*
	if (deltaS<0) return true;
	else{
		double coin = xx.rand();

		if (coin<exp(-deltaS)) return true;
		else return false;
	}
	*/
}


double Metropolis_for_one_link(int site, int dir, double tuner, dc* ufield_new) {
	// This is an ULTRA LOCAL Metropolis change on the gauge variables, given one link
	// Generate a random SU(2) matrix X=1+tuner*R,
	// Generate proposal as U(mu,x) ---> X.U(mu,x)
	
	// generate SU(2) random matrix X=1-tuner*R,
	dc* X = new dc[4];
	for(int i=0; i<Ncolsquare; i++) X[i] = tuner*dc(.5-xx.rand(),.5-xx.rand());
	for(int i=0; i<Ncolsquare; i+=Ncol_plus_one) X[i] += dc(1.,0.);
	norm_su2(X);

	// Generate proposal
	dc Uold[4], Unew[4];
	for(int i=0; i<4; i++) Uold[i] = ufield[(dir*nsites+site)*Ncolsquare+i];
	mult_C_equals_AB_for_SU2(Unew,X,Uold);

	// CALCULATE STAPLE
	dc bin1[4], bin2[4], bin3[4], staplesum[4];
	for(int i=0; i<4; i++) staplesum[i] = dc(0.,0.);

	for(int nu=0; nu<dim; nu++) {
		if(nu!=dir) {
			// pstaple(bin1,site,dir,nu,0);
			// nstaple(bin1,site,dir,nu,0);
			// for(int i=0; i<4; i++) staplesum[i] += bin1[i];

			// Calculate UPPER staple
			dc up[4], left[4], down[4];

			int n_plus_mu = neighbor_plus[dir*nsites+site];
			int n_plus_nu = neighbor_plus[nu*nsites+site];
			
			for(int i=0; i<4; i++) up[i]   = ufield[(nu*nsites+n_plus_mu)*Ncolsquare+i];  // U_nu(n+mu)
			for(int i=0; i<4; i++) left[i] = ufield[(dir*nsites+n_plus_nu)*Ncolsquare+i];  // U_mu(n+nu)
			for(int i=0; i<4; i++) down[i] = ufield[(nu*nsites+site)*Ncolsquare+i];       // U_nu(n)

			mult_C_equals_ABdagger_for_SU2(bin1,up,left);
			mult_C_equals_ABdagger_for_SU2(bin2,bin1,down);


			// Calculate UPPER staple
			int n_minus_nu = neighbor_minus[nu*nsites+site];
			int n_minus_nu_plus_mu = neighbor_plus[dir*nsites+n_minus_nu];
			
			for(int i=0; i<4; i++) down[i] = ufield[(nu*nsites+n_minus_nu_plus_mu)*Ncolsquare+i];  // U_nu(n+mu-nu)
			for(int i=0; i<4; i++) left[i] = ufield[(dir*nsites+n_minus_nu)*Ncolsquare+i];          // U_mu(n-nu)
			for(int i=0; i<4; i++) up[i]   = ufield[(nu*nsites+n_minus_nu)*Ncolsquare+i];          // U_nu(n-mu)

			mult_C_equals_AdaggerB_for_SU2(bin1,left,up);
			mult_C_equals_AdaggerB_for_SU2(bin3,down,bin1);


			for(int i=0; i<4; i++) staplesum[i] += bin2[i] + bin3[i];
		}
	}

	// Calculate deltaS
	for(int i=0; i<4; i++) bin1[i] = Unew[i] - Uold[i];
	mult_C_equals_AB_for_SU2(bin2,bin1,staplesum);

	double deltaS = -beta/2.*real(bin2[0]+bin2[3]);

	// Update exiting proposal
	for(int i=0; i<4; i++) ufield_new[(dir*nsites+site)*Ncolsquare+i] = Unew[i];

	return deltaS;
}


// This performs a sweep (== update every link variable once) and then perform metropolis
bool Metropolis_sweep_gauge(double tuner) {
	// Create a bin and fill with gauge variables
	
	dc *proposal = new dc[ufielddimension];
	for(int ii=0; ii<ufielddimension; ii++) proposal[ii] = ufield[ii];
	
	// Perform the sweep
	double deltaS = 0.;

	for(int site=0; site<nsites; site++)
	for(int dir=0; dir<dim; dir++) {
		// int site = xx.randInt(nsites);
		// int dir = xx.rand(dim);
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

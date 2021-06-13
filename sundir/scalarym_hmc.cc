// THIS IS TO USE HMC FOR GAUGE FIELDS
double get_Hamiltonian(double* P) {
    double H = get_wilson_action();

    for(int i=0; i<(nsites*dim*(Ncolsquare-1)); i++) H += .5*P[i]*P[i]; 

    return H;
}


// This generates gaussian distributed momenta
void refresh_mom(double * mom_comp ){
    int dim_mom_comp = (Ncolsquare-1)*nsites*dim;
    for(int i=0; i<dim_mom_comp; i++ ){ 
        mom_comp[i] = xx.randNorm(0,1);
    }
}

void compute_forces( double* force_comp){    
    dc aux_pauli[Ncolsquare], aux1[Ncolsquare], U[Ncolsquare], staplesum[Ncolsquare],  aux_staple[Ncolsquare];
    //Initialize bins
    
    for(int n=0; n<nsites; n++)
    for(int mu=0; mu<dim; mu++)
    for(int a=0; a<(Ncolsquare-1); a++){
        // set staplesum to zero
        for(int i=0; i<Ncolsquare; i++) staplesum[i] = dc(0.,0.);
        
        // get sum staple
        for(int nu=0; nu<dim; nu++){
            if(nu!=mu){
                pstaple(aux_staple, n, mu, nu, 0);
			    nstaple(aux_staple, n, mu, nu, 0);
			    for(int i=0; i<Ncolsquare; i++) staplesum[i] += aux_staple[i];
            }
        }

        // get pauli and gauge fields of interest
        for(int i=0; i<Ncolsquare; i++) aux_pauli[i] = .5*pauli[a*Ncolsquare + i]; 
        for(int i=0; i<Ncolsquare; i++) U[i] = ufield[(mu*nsites+n)*Ncolsquare +i]; 

        // multiply pauli*gauge*sumstaple
        mult_C_equals_AB_for_SU2(aux1, aux_pauli, U); 
        mult_C_equals_AB_for_SU2(aux_staple, aux1, staplesum); 
        
        force_comp[(mu*nsites + n)*(Ncolsquare-1)+a] = beta / float(Ncol) * imag(aux_staple[0] + aux_staple[3]);
    }

    return;
}


// This calculates I_U(epsilon)
void leap_U(double epsilon, dc* U, double* P) {
    double C = cos(epsilon);
    double S = sin(epsilon);

    for(int site=0; site<nsites; site++)
    for(int mu=0; mu<dim; mu++) {
        
        // Calculate exp(i*epsilon*P_mu(n)) giving the array of components of gauge momenta {P_\mu^a(n)}
        double P1 = P[(mu*nsites+site)*3];
        double P2 = P[(mu*nsites+site)*3+1];
        double P3 = P[(mu*nsites+site)*3+2];

        dc expP[4];
        expP[0] = dc(C,S*P3);
        expP[1] = dc(P2,P1)*S;
        expP[2] = dc(-P2,P1)*S;
        expP[3] = dc(C,-S*P3);

        // Produce candidates U' = expP*U
        dc uold[4], uprime[4];
        for (int i=0; i<Ncolsquare; i++) uold[i] = U[(mu*nsites+site)*Ncolsquare+i];
        mult_C_equals_AB_for_SU2(uprime,expP,uold);

        // Update leaped U
        for (int i=0; i<Ncolsquare; i++) U[(mu*nsites+site)*Ncolsquare+i] = uprime[i];
    }

    return;
}

void leap_P(double epsilon, double* P) {
    double* forces = new double[(Ncolsquare-1)*nsites*dim];
    compute_forces(forces);

    for(int site=0; site<nsites; site++)
    for(int mu=0; mu<dim; mu++)
    for(int a=0; a<(Ncolsquare-1); a++) {
        P[(mu*nsites+site)*(Ncolsquare-1)+a] -= epsilon*forces[(mu*nsites+site)*(Ncolsquare-1)+a];
    }

    delete[] forces;
    return;
}



bool jump_HMC(double epsilon, double tau) {

    int Njump = int(tau/epsilon);

    // Copy field configurations and initialize momenta
    dc* ufield_old = new dc[ufielddimension];
    for(int i=0; i<ufielddimension; i++) ufield_old[i] = ufield[i];

    double* momenta = new double[nsites*dim*(Ncolsquare-1)];
    refresh_mom(momenta);


    double deltaH = -get_Hamiltonian(momenta);


    for(int jump=0; jump<Njump; jump++) {
        leap_P(.5*epsilon,momenta);
        leap_U(epsilon,ufield,momenta);
        leap_P(.5*epsilon,momenta);
    }


    // METROPOLIS COINFLIP
    deltaH += get_Hamiltonian(momenta);

    bool coin = boltzmann_coin_flip(deltaH);
    delete[] momenta;


    if (coin==false) {
        for(int i=0; i<ufielddimension; i++) ufield[i] = ufield_old[i];
        cout << "Reject\n";
        return false;
    }
    else {
        cout << "         Accept\n";
        return true;
    }

}
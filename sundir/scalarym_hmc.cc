
 

    



// STEP 1 : GENERATE MOMENTA

void refresh_mom(double * mom_comp ){
    int dim_mom_comp = (Ncolsquare-1)*nsites*dim;
    for(int i=0; i<dim_mom_comp; i++ ){ 
        mom_comp[i] = xx.randNorm(0,1);
    }
}

//double* force_comp = new double[];

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
        for(int i=0; i<Ncolsquare; i++) aux_pauli[i] = pauli[a*Ncolsquare + i]; 
        for(int i=0; i<Ncolsquare; i++) U[i] = ufield[(mu*nsites+n)*Ncolsquare +i]; 

        // multiply pauli*gauge*sumstaple
        mult_C_equals_AB_for_SU2(aux1, aux_pauli, U); 
        mult_C_equals_AB_for_SU2(aux_staple, aux1, staplesum); 
        
        force_comp[(mu*nsites + n)*(Ncolsquare-1)+a] = beta / float(Ncol) * imag(aux_staple[0] + aux_staple[3]);

    }

}

/*
void 


int force_dim = nsites * dim * (Ncolsquare-1);
dc *ym_forces = new dc[force_dim];


void compute_forces(int site, int dir, dc*forces){

    for(int i=0; i<force_dim; i++){
        forces[i] = 
    }
    

}
*/
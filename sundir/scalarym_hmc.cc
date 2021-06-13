
 /*  
    dc *sigma1 = new dc[4];
    sigma1[1] = dc(0.,0.);
    sigma1[2] = dc(1.,0.);
    sigma1[3] = dc(1.,0.);
    sigma1[4] = dc(0.,0.);
    
    dc *sigma2 = new dc[4];
    sigma2[1] = dc(0.,0.);
    sigma2[2] = dc(0.,-1.);
    sigma2[3] = dc(0.,1.);
    sigma2[4] = dc(0.,0.);
    
    dc *sigma3 = new dc[4];
    sigma3[1] = dc(1.,0.);
    sigma3[2] = dc(0.,0.);
    sigma3[3] = dc(0.,0.);
    sigma3[4] = dc(-1.,0.);
    */



// STEP 1 : GENERATE MOMENTA

void refresh_mom(double * mom_comp ){
    int dim_mom_comp = (Ncolsquare-1)*nsites*dim;
    for(int i=0; i<dim_mom_comp; i++ ){ 
        mom_comp[i] = xx.randNorm(0,1);
    }
}

void 

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
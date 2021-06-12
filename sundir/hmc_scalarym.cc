{
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
}
int get_site(){

}


int force_dim = nsites * dim * (Ncolsquare-1);
dc *ym_forces = new dc[force_dim];


void compute_forces(int site, int dir, dc*forces){

    for(int i=0; i<force_dim; i++){
        forces[i] = 
    }
    

}
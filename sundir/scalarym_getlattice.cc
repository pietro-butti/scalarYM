int get_site_index(int t, int x, int y){
    return y + ny * (x+nx * t);
}

int get_link_index(int site, int mu ){
    return mu*nsites + site;
}

/* Given a link index, get_ufield stores in *res
the value of the associated gauge field */
void get_ufield(int link, dc* res){
    for(int i=0; i<Ncolsquare; i++ ) res[i] = ufield[link*Ncolsquare +i];
}

/* Given a site and two directions, 
this function stores in *res the plaquette indeces*/
void get_plaq_index(int site, int mu, int nu, int* res){
    res[0] = get_link_index(site, mu);
    res[1] = get_link_index(neighbor_plus[mu*nsites + site], nu);
    res[2] = get_link_index(neighbor_plus[nu*nsites + site], mu);
    res[3] = get_link_index(site, nu);
}


double get_wilson_action(){

    double action = 0.;    
    for(int site=0; site<nsites; site++)
    for(int nu=0; nu<dim; nu++)
<<<<<<< HEAD
    for(int mu=0; mu<nu; mu++)if(mu!=nu){

=======
    for(int mu=0; mu<nu; mu++) if(mu!=nu) 
    {
>>>>>>> origin/master
        int* plaq_index = new int[4];
        get_plaq_index(site, mu, nu, plaq_index);
        
        dc * aux1 = new dc[4];
        dc * aux2 = new dc[4];
        
        dc* a = new dc[4]; get_ufield( plaq_index[0],a);
        dc* b = new dc[4]; get_ufield( plaq_index[1],b);
        dc* c = new dc[4]; get_ufield( plaq_index[2],c);
        dc* d = new dc[4]; get_ufield( plaq_index[3],d);


        mult_C_equals_AB_for_SU2(aux1, a, b);
        mult_C_equals_ABdagger_for_SU2(aux2, aux1, c);
        mult_C_equals_ABdagger_for_SU2(aux1, aux2, d);

<<<<<<< HEAD
        //for(int i=0; i<4; i++) cout << aux1[i]<<endl;       
        action += 1. - 1./2. * real(aux1[0] + aux1[3]); 
        // action += 1./double(Ncol*nsites*dim*(dim-1))*real(aux1[0] + aux1[3]);
    }
    return action*beta;
=======
        action += 1. - 1./2. * real(aux1[0] + aux1[3]); 
        // action += 2./double(Ncol*nsites*dim*(dim-1))*real(aux1[0] + aux1[3]);
    }
    return 2.*action*beta;
>>>>>>> origin/master
}


double get_plaquette(){

    dc bin1[4], bin2[4];
    double S = 0.;

    for(int site=0; site<nsites; site++)
    for(int nu=0; nu<dim; nu++)
    for(int mu=0; mu<nu; mu++) {
        
        dc right[4]; dc up[4]; dc left[4]; dc down[4];
        for(int i=0; i<Ncolsquare; i++) right[i] = ufield[(mu*nsites+site)*Ncolsquare+i];
        for(int i=0; i<Ncolsquare; i++) up[i]    = ufield[(nu*nsites+neighbor_plus[mu*nsites+site])*Ncolsquare+i];
        for(int i=0; i<Ncolsquare; i++) left[i]  = ufield[(mu*nsites+neighbor_minus[nu*nsites+site])*Ncolsquare+i];
        for(int i=0; i<Ncolsquare; i++) down[i]  = ufield[(nu*nsites+site)*Ncolsquare+i];

        mult_C_equals_AB_for_SU2(bin1,right,up);
        mult_C_equals_ABdagger_for_SU2(bin2,bin1,left);
        mult_C_equals_ABdagger_for_SU2(bin1,bin2,down);

        S += 1./double(Ncol*nsites*dim*(dim-1))*real(bin1[0]+bin1[3]);        
    }
    
    return  S;
}


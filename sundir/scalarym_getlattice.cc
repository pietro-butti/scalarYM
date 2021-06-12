int get_site_index(int t, int x, int y){
    return y + ny * (x+nx * t);
}

int get_link_index(int site, int mu ){
    return mu*nsites + site;
}

void get_ufield(int link, dc* res){
    for(int i=0; i<Ncolsquare; i++ ) res[i] = ufield[link*Ncolsquare +i];
}

void get_plaq_index(int site, int mu, int nu, int* res){
    res[0] = get_link_index(site, mu);
    res[1] = get_link_index(neighbor_plus[mu*nsites + site], nu);
    res[2] = get_link_index(neighbor_plus[nu*nsites + site], mu);
    res[3] = get_link_index(site, nu);
}

double get_wilson_action(){

    double action = 0.;
    for(int site=0; site<nsites; site++)
    for(int nu=0; nu<4; nu++)
    for(int mu=0; mu<nu; mu++){
        
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

        //action += 1. - 1./2. * real(aux1[0] + aux1[3]); 
        action += real(aux1[0] + aux1[3]);
    }
    return action/(nsites*6.);
}
#if Ncol==7

// Representation 21:
temporary_trace[1]=********weyl
zero_mom_loops_now[1][x]+=temporary_trace[1];

// Representation 28:
temporary_trace[2]=f*f-temporary_trace[1];
zero_mom_loops_now[2][x]+=temporary_trace[2];

// Representation 35:
temporary_trace[3]=********weyl
zero_mom_loops_now[3][x]+=temporary_trace[3];

// Representation 48:
temporary_trace[4]=f*fbar-1.
zero_mom_loops_now[4][x]+=temporary_trace[4];

/**********************************************************/
/* Now we calculate the 112, to use it for the 84: */
      aux=f*temporary_trace[1]-temporary_trace[3];
/**********************************************************/
    
// Representation 84:
temporary_trace[5]=temporary_trace[2]*f-aux;
zero_mom_loops_now[5][x]+=temporary_trace[5];

// Representation 112:
temporary_trace[6]=aux;
zero_mom_loops_now[6][x]+=temporary_trace[6];

// Representation 140:
temporary_trace[7]=conj(temporary_trace[1])*f-fbar;
zero_mom_loops_now[7][x]+=temporary_trace[7];

// Representation 189:
temporary_trace[8]=temporary_trace[4]*f-conj(temporary_trace[7])-fbar;
zero_mom_loops_now[8][x]+=temporary_trace[8];

/**********************************************************/
/* Now we calculate the 210 and the 210', to use them for the 196: */
      aux=temporary_trace[3]*f-conj(temporary_trace[3]);
      aux2=****weyl;
/**********************************************************/
    
// Representation 196:
temporary_trace[9]=f*(temporary_trace[6]-temporary_trace[5])-aux+aux2;
zero_mom_loops_now[9][x]+=temporary_trace[9];

// Representation 210:
temporary_trace[10]=aux;
zero_mom_loops_now[10][x]+=temporary_trace[10];

// Representation 210':
temporary_trace[11]=aux2;
zero_mom_loops_now[11][x]+=temporary_trace[11];



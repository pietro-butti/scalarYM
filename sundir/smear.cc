void smear() {
  
  int site, dir, i, iteration, index;
  int first_entry, second_entry, projection;
  double v0, v1, v2, v3, inverse_vmodulus;
  double real_trace_variation, real_trace_before;
  
  for (i=0; i<ufielddimension; i++) {
    smeared_ufield[i]=ufield[i];
  }
  
  for (iteration=0; iteration<smearing_iterations; iteration++) {
    for (site=0; site<nsites; site++) 
    for (dir=first_smearing_direction; dir<dim; dir++) {
      pstaple( staple, site, dir, first_smearing_direction, 1);
      nstaple( staple, site, dir, first_smearing_direction, 1);
      for (i=0; i<Ncolsquare; i++) {
        u1[i]=link_smearing_coeff*smeared_ufield[(dir*nsites+site)*Ncolsquare+i]+staples_smearing_coeff*staple[i];
// u1 is initially defined as the linear combination 
// of the link and the staples
        u3[i]=dc(0.,0.);
      }
      for (i=0; i<Ncolsquare; i+=Ncol_plus_one) {
        u3[i]=dc(1.,0.);
      }
      
      real_trace_before=-1000.;
      real_trace_variation=real_trace_before;
      projection=0;
      
      while ((projection<max_number_of_projections) && (fabs(real_trace_variation)>convergence_tolerance)) {
        projection++;
        for (first_entry=0;first_entry<Ncol-1;first_entry++)
        for (second_entry=first_entry+1;second_entry<Ncol;second_entry++) {
          for (i=0; i<Ncolsquare; i++) {
            tempor2[i]=dc(0.,0.);
          }
          for (i=0; i<Ncolsquare; i+=Ncol_plus_one) {
            tempor2[i]=dc(1.,0.);
          }
          
          v0 = real( u1[first_entry*Ncol+first_entry]
                   + u1[second_entry*Ncol+second_entry] );
          v1 = imag( u1[first_entry*Ncol+second_entry]
                   + u1[second_entry*Ncol+first_entry] );
          v2 = real( u1[first_entry*Ncol+second_entry]
                   - u1[second_entry*Ncol+first_entry] );
          v3 = imag( u1[first_entry*Ncol+first_entry]
                   - u1[second_entry*Ncol+second_entry] );
// 
// The v_mu's are (up to an irrelavant factor 1/2, which is
// reabsorbed by the normalization) the real parts of the 
// coefficients in the decomposition of the (2x2) block of
// u1 as:
// 
//   u1 = v0 + i \sigma_k v_k
// 
// Note that:
// 1) since U is in SU(N), Re Tr (U V\dagger) depends only on 
// the real part of the v_mu's;
// 2) Re Tr (U V\dagger) is maximized by the same U which also
// maximizes Re Tr (U B\dagger), where the coefficients of
// B are real and normalized to 1 (so that B is in SU(2));
// 3) Re Tr (U B\dagger) is maximized by U=B.
// 
          inverse_vmodulus=1./sqrt(v0*v0+v1*v1+v2*v2+v3*v3);
          v0*=inverse_vmodulus;
          v1*=inverse_vmodulus;
          v2*=inverse_vmodulus;
          v3*=inverse_vmodulus;
          tempor2[first_entry*Ncol+first_entry]=dc(v0,v3);
          tempor2[first_entry*Ncol+second_entry]=dc(v2,v1);
          tempor2[second_entry*Ncol+first_entry]=dc(-v2,v1);
          tempor2[second_entry*Ncol+second_entry]=dc(v0,-v3);
          mult_C_equals_ABdagger(tempor,tempor2,u1);
// We calculate the new value of the trace of Re Tr (U V\dagger):
  
          v0=0.;
          for (i=0; i<Ncolsquare; i+=Ncol_plus_one) {
            v0+=real(tempor[i]);
          }
//   cout << "ReTr (U Vdag) = " << v0 << endl;
          real_trace_variation=v0-real_trace_before;
          if (real_trace_variation<-convergence_tolerance) {
            cout << "# The trace is decreasing!\n";
            printf("# real_trace_variation=%e\n",real_trace_variation);
//             do {;} while (1);
          }
          real_trace_before=v0;
	  
// Finally, we prepare for the next SU(2) subgroup, by modifying
// u1 to keep track of the SU(2) matrix that we just found
// - note that, at the next iteration, we want to maximize 
// Re Tr (U_2 U_1 V\dagger), which can be rewritten as
// Re Tr (U_2 W\dagger), if we set W = V U_1\dagger -, 
// and by multiplying U_1 by u3 (after which, we
// redefine u3 to be the product)
// 
          mult_C_equals_ABdagger(tempor,u1,tempor2);
          mult_C_equals_AB(u2,tempor2,u3);
          for (i=0; i<Ncolsquare; i++) {
            u1[i]=tempor[i];
            u3[i]=u2[i];
          }
          
        }
// End of the cycle over SU(2) subgroups        
        
      }
// End of the trace maximization cycle
      
      if (projection>=max_number_of_projections) {
        cout << "# Failed to converge\n";
        cout << "#  |real_trace_variation|="<<fabs(real_trace_variation)<<endl;
//             do {;} while (1);
      }
//       else {
//         cout << "Converged in " << projection << " steps" << endl;
//       }
            
      mult_C_equals_AB(u2,tempor2,u3);
      
      for (i=0; i<Ncolsquare; i++) {
        tempfield[(dir*nsites+site)*Ncolsquare+i]=u2[i];
      }
      
    }
// End of the cycle over the links
    
    for (site=0; site<nsites; site++) 
    for (dir=first_smearing_direction; dir<dim; dir++) {
      index=(dir*nsites+site)*Ncolsquare;
      for (i=0; i<Ncolsquare; i++) {
        smeared_ufield[index+i]=tempfield[index+i];
      }
    }
    
    
// Unitarity check: 
// 
//     for (site=0; site<nsites; site++) 
//     for (dir=0; dir<dim; dir++) {  
//       for (i=0; i<n; i++)
//       for (j=0; j<n; j++) {
//         u3[i][j]=smeared_ufield[i][j][site][dir];
//         u4[i][j]=smeared_ufield[i][j][site][dir];
//       }
//       mult_C_equals_ABdagger(tempor,u3,u4);
//       print_matrix(tempor);
//       mult_C_equals_AdaggerB(tempor,u3,u4);
//       print_matrix(tempor);
//       if (n==3) {
//         dc determinante=
//             u3[0][0]*u3[1][1]*u3[2][2]
//           + u3[0][1]*u3[1][2]*u3[2][0]
//           + u3[0][2]*u3[1][0]*u3[2][1]
//           - u3[2][0]*u3[1][1]*u3[0][2]
//           - u3[2][1]*u3[1][2]*u3[0][0]
//           - u3[2][2]*u3[1][0]*u3[0][1];
//         printf("iteration=%d site=%d dir=%d  det=(%12.10lf +I*(%12.10lf))\n",
//           iteration, site, dir,
//           real(determinante),
//           imag(determinante));
//       }
//     }
// 

  }
// End of the cycle over the smearing iterations

}

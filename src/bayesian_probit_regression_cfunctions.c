#include <math.h>
#include <stdlib.h>

void RandNorm(int *num_samples, int *num_unif_samples, double *mean,
                   double *variance, double *samples){


  double std_dev = sqrt(*variance);
  int nsamp = *num_samples;
  int nunifsamp = *num_unif_samples;

  double sqrt_b = sqrt(12.0 * (double)nunifsamp);

  int i, j;

  for(i=0; i<nsamp; i++){
    double unif_sum = 0.0;
    for(j=0; j<nunifsamp; j++){
      double curr_rand = (double)rand() / (double)RAND_MAX;
      unif_sum += curr_rand * sqrt_b;
    }

    double curr_mean_std = (unif_sum / (double)nunifsamp) - (sqrt_b / 2.0);

    double curr_sample = curr_mean_std * std_dev + (*mean);

    *(samples + i) = curr_sample;

  }

}

void RandTruncNorm(int *num_samples, int *num_unif_samples, double *mean, double *variance,
                   double *cutoff, int *greater, double *samples, int *attempts_failed){


  int i = 0;

  while(i < *num_samples){

    int num_attempts = 0;

    double curr_sample[1];
    int pointer_to_one[1] = {1};

    while(num_attempts < 100000){

      RandNorm(pointer_to_one, num_unif_samples, mean, variance, curr_sample);
      if(*greater == 1){
        num_attempts++;
        if(*curr_sample > *cutoff){
          *(samples + i) = *curr_sample;
          goto next_iter;
        }
      }
      else{
        num_attempts++;
        if(*curr_sample <= *cutoff){
          *(samples + i) = *curr_sample;
          goto next_iter;
        }
      }
    }

    if(num_attempts >= 100000){
      *attempts_failed = 1;
      return;
    }

    next_iter:

    i++;
  }

}



void CreateDiagonalMatrix(double *A, int *num_rows){
  int i, j;
  int nrow = *num_rows;
  for(i=0;i<nrow;i++){
    for(j=0;j<nrow;j++){
      if(i==j){
        *(A+i+j*nrow) = 1.0;
      }
      else{
        *(A+i+j*nrow) = 0.0;
      }
    }
  }
}



void CreateZeroMatrix(double *A, int *num_rows, int *num_cols){
  long i, j;
  int nrow = *num_rows;
  int ncol = *num_cols;
  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      *(A+i+j*nrow) = 0.0;
    }
  }
}

void CopyMatrixToAnother(double *A, double *B, int *num_rows, int *num_cols){
  long i, j;
  int nrow = *num_rows;
  int ncol = *num_cols;
  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      *(B+i+j*nrow) = *(A+i+j*nrow);
    }
  }
}

void TransposeMatrix(double *A, int *num_rows, int *num_cols, double *tA){
  long i, j;
  int nrow = *num_rows;
  int ncol = *num_cols;
  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      *(tA+j+i*ncol) = *(A+i+j*nrow);
    }
  }
}


void AddOrSubtractMatrices(double *A, double *B, int *num_rows,
                           int *num_cols,  double *sum, int *is_addition){
  int i, j;
  int nrow = *num_rows;
  int ncol = *num_cols;
  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      if(*is_addition == 1){
        *(sum+i+j*nrow) = (*(A+i+j*nrow)) + (*(B+i+j*nrow));
      }
      else{
        *(sum+i+j*nrow) = (*(A+i+j*nrow)) - (*(B+i+j*nrow));
      }

    }
  }
}


void MultiplyMatrices(double *A, double *B, int *num_rows_A,
                      int *num_cols_A,  int *num_cols_B, double *AB){

  int i, j;
  int nrow_A = *num_rows_A;
  int ncol_A = *num_cols_A;
  int ncol_B = *num_cols_B;
  int nrow_B = ncol_A;
  for(i=0;i<nrow_A;i++){
    for(j=0;j<ncol_B;j++){
      int k;
      double curr_row_col_prod = 0.0;
      for(k=0;k<ncol_A;k++){
        curr_row_col_prod += (*(A+i+k*nrow_A)) * (*(B+j*nrow_B+k));
      }
      *(AB+i+j*nrow_A) = curr_row_col_prod;
    }
  }
}


void MultiplyMatrixWithConstant(double *A, double *c, int *num_rows, int *num_cols, double *cA){
  int i, j;
  int nrow = *num_rows;
  int ncol = *num_cols;
  for(i=0;i<nrow;i++){
    for(j=0;j<ncol;j++){
      *(cA+i+j*nrow) = (*c) * (*(A + i + j*nrow));
    }
  }
}



void GaussEliminationInverse(double *A, int *num_rows, double *Ainv, int *inverse_failed){

  int nrow = *num_rows;
  int ncol = nrow;
  double copyA[nrow*nrow];
  int i, j;

  CopyMatrixToAnother(A, copyA, num_rows, num_rows);
  CreateDiagonalMatrix(Ainv, num_rows);

  for(i=0;i<nrow;i++){
    int k;
    double curr_pivot = (*(copyA+i+i*nrow));

    if(curr_pivot == 0){
      *inverse_failed = 1;
      return;
    }

    for(k=0;k<ncol;k++){
      *(copyA+i+k*nrow) /= curr_pivot;
      *(Ainv+i+k*nrow) /= curr_pivot;
    }

    for(j=i+1;j<nrow;j++){
      int l;
      double curr_eliminator = (*(copyA+j+i*nrow));
      for(l=0;l<ncol;l++){
        *(copyA+j+l*nrow) -= curr_eliminator * (*(copyA+i+l*nrow));
        *(Ainv+j+l*nrow) -= curr_eliminator * (*(Ainv+i+l*nrow));
      }
    }
  }

  for(i=(nrow-1);i>0;i--){
    for(j=i-1;j>=0;j--){
      int l;
      double curr_backward_eliminator = (*(copyA+j+i*nrow));
      for(l=0;l<ncol;l++){
        *(copyA+j+l*nrow) -= curr_backward_eliminator * (*(copyA+i+l*nrow));
        *(Ainv+j+l*nrow) -= curr_backward_eliminator * (*(Ainv+i+l*nrow));
      }
    }
  }
}


void CholeskyDecomposition(double *A, double *L, int *num_rows,
                           int *positive_deifinte_failed){
  int nrow = *num_rows;
  int i, j;
  int pointer_to_one[1] = {1};

  double copyA[nrow*nrow];
  CopyMatrixToAnother(A, copyA, num_rows, num_rows);

  double copyL[nrow*nrow];

  CreateDiagonalMatrix(L, num_rows);

  for(i=0; i < nrow; i++){
    CreateDiagonalMatrix(copyL, num_rows);

    double curr_pivot = *(copyA + i + i*nrow);
    if(curr_pivot <= 0){
      *positive_deifinte_failed = 1;
      return;
    }

    double sqrt_pivot = sqrt(curr_pivot);

    for(j=i; j < nrow; j++){
      *(copyL + j + i*nrow) = (*(copyA + j + i*nrow))/sqrt_pivot;
    }

    *(copyA + i + i * nrow) = 1.0;
    if(i < (nrow - 1)){

      double curr_column[nrow];
      CreateZeroMatrix(curr_column, num_rows, pointer_to_one);

      for(j=(i+1); j < nrow; j++){
        *(curr_column + j) = *(copyA + j + i*nrow);
      }

      for(j=(i+1); j < nrow; j++){
        *(copyA + i + j*nrow) = 0.0;
        *(copyA + j + i*nrow) = 0.0;
      }

      double cross_prod_column[nrow*nrow];
      MultiplyMatrices(curr_column, curr_column, num_rows, pointer_to_one,
                       num_rows, cross_prod_column);

      double c[1] = {1/curr_pivot};
      double normalized_cross_prod[nrow*nrow];
      MultiplyMatrixWithConstant(cross_prod_column, c, num_rows,
                                 num_rows, normalized_cross_prod);

      double copyAtemp[nrow * nrow];
      int is_addition[1] = {0};
      AddOrSubtractMatrices(copyA, normalized_cross_prod, num_rows, num_rows,
                            copyAtemp, is_addition);

      CopyMatrixToAnother(copyAtemp, copyA, num_rows, num_rows);

    }

    double oldL[nrow * nrow];
    CopyMatrixToAnother(L, oldL, num_rows, num_rows);
    MultiplyMatrices(oldL, copyL, num_rows, num_rows, num_rows, L);

  }
}


void BayesProbitC(int *y, double *X, double *covar_prior, int *num_variables,
                  int *num_observations, int *num_mcmc, double *coeff_samples,
                  int *positive_deifinte_failed, int *attemps_failed, int *inverse_failed){

  int nvar = *num_variables;
  int nobs = *num_observations;
  int nchain = *num_mcmc;

  int i,j,k;
  int pointer_to_one[1] = {1};
  int pointer_to_nunif[1] = {50};
  double cut_off[1] = {0.0};

  double X_transpose[nvar*nobs];
  TransposeMatrix(X, num_observations, num_variables, X_transpose);

  double Xcrossprod[nvar*nvar];
  MultiplyMatrices(X_transpose, X, num_variables, num_observations,
                   num_variables, Xcrossprod);

  double prec_prior[nvar*nvar];
  GaussEliminationInverse(covar_prior, num_variables, prec_prior, inverse_failed);

  if(*inverse_failed == 1){
    return;
  }

  double Vinv[nvar*nvar];
  int is_addition[1] = {1};
  AddOrSubtractMatrices(Xcrossprod, prec_prior, num_variables, num_variables,
                        Vinv, is_addition);

  double V[nvar*nvar];
  GaussEliminationInverse(Vinv, num_variables, V, inverse_failed);

  if(*inverse_failed == 1){
    return;
  }

  double L[nvar*nvar];

  CholeskyDecomposition(V, L, num_variables, positive_deifinte_failed);

  if(*positive_deifinte_failed == 1){
    return;
  }

  double S[nvar*nobs];
  MultiplyMatrices(V, X_transpose, num_variables, num_variables, num_observations, S);

  double H[nobs];
  double W[nobs];
  double Q[nobs];

  for(j=0; j<nobs; j++){

    double Xj[nvar];
    double Sj[nvar];

    for(k=0; k<nvar; k++){
      *(Xj + k) = *(X + j + nobs*k);
      *(Sj + k) = *(S + k + nvar*j);
    }

    double Hj[1];
    MultiplyMatrices(Xj, Sj, pointer_to_one, num_variables, pointer_to_one, Hj);
    *(H + j) = *(Hj);

    double Wj[1];
    *(Wj) = (*Hj) / (1 - (*Hj));
    *(W + j) = *(Wj);

    double Qj[1];
    *(Qj) = *(Wj) + 1.0;
    *(Q + j) = *(Qj);

  }

  double Z[nobs];

  for(i=0; i<nobs; i++){
    double Zi[1];
    double mean_zi[1] = {0.0};
    double var_zi[1] = {1.0};

    if((*(y + i)) == 0){
      int greater[1] = {0};
      RandTruncNorm(pointer_to_one, pointer_to_nunif,
                    mean_zi, var_zi, cut_off, greater, Zi, attemps_failed);
    }
    else{
      int greater[1] = {1};
      RandTruncNorm(pointer_to_one, pointer_to_nunif,
                    mean_zi, var_zi, cut_off, greater, Zi, attemps_failed);
    }

    if(*attemps_failed == 1){
      return;
    }

    *(Z + i) = *(Zi);
  }

  double B[nvar];
  MultiplyMatrices(S, Z, num_variables, num_observations, pointer_to_one, B);

  for(i = 0; i<nchain; i++){
    for(j = 0; j<nobs; j++){

      double z_old;
      z_old = *(Z + j);

      double mj[1];
      double mj_new[1];
      double Xj[nvar];
      double Sj[nvar];

      for(k=0; k<nvar; k++){
        *(Xj + k) = *(X + j + nobs*k);
        *(Sj + k) = *(S + k + nvar*j);
      }

      MultiplyMatrices(Xj, B, pointer_to_one, num_variables, pointer_to_one, mj);

      *(mj_new) = *(mj) - (*(W + j)) * ((*(Z + j)) - (*mj));

      double zj[1];
      double qj[1];

      *(qj) = *(Q + j);

      if((*(y + j)) == 0){
        int greater[1] = {0};
        RandTruncNorm(pointer_to_one, pointer_to_nunif,
                      mj_new, qj, cut_off, greater, zj, attemps_failed);
      }
      else{
        int greater[1] = {1};
        RandTruncNorm(pointer_to_one, pointer_to_nunif,
                      mj, qj, cut_off, greater, zj, attemps_failed);
      }

      if(*attemps_failed == 1){
        return;
      }

      *(Z + j) = *(zj);

      double Bnew[nvar];
      double z_new[1];
      double zS[nvar];

      *(z_new) = *(zj) - z_old;
      MultiplyMatrixWithConstant(Sj, z_new, num_variables, pointer_to_one, zS);

      AddOrSubtractMatrices(B, zS, num_variables, pointer_to_one, Bnew, is_addition);

      CopyMatrixToAnother(Bnew, B, num_variables, pointer_to_one);

    }

    double T[nvar];

    for(k=0; k<nvar; k++){
      double tk[1];
      double mean_k[1] = {0.0};
      double var_k[1] = {1.0};

      RandNorm(pointer_to_one, pointer_to_nunif, mean_k, var_k, tk);

      *(T + k) = *(tk);
    }

    double LT[nvar];
    double BLT[nvar];
    MultiplyMatrices(L, T, num_variables, num_variables, pointer_to_one, LT);
    AddOrSubtractMatrices(B, LT, num_variables, pointer_to_one, BLT, is_addition);

    for(k=0; k<nvar; k++){
      *(coeff_samples + k + nvar*i) = *(BLT + k);
    }

  }

}

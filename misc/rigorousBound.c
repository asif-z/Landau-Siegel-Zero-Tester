#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <flint/arb.h>
#include <time.h>
#include<unistd.h>

// sets output = |a + b*e^{i\theta}|^2
void norm(arb_t output, arb_t a, arb_t b, arb_t theta, long prec) {
    arb_t real;
    arb_t imag;
    arb_init(real);
    arb_init(imag);

    arb_cos(real, theta, prec);
    arb_mul(real, real, b, prec);
    arb_add(real, real, a, prec);
    arb_mul(real, real, real, prec);

    arb_sin(imag, theta, prec);
    arb_mul(imag, imag, b, prec);
    arb_mul(imag, imag, imag, prec);

    arb_add(output, real, imag, prec);
}

int main(int argc, char** argv)
{
    long prec = 70;
    arb_t sigma;
    arb_t logQ;
    arb_t r;
    arb_t lambda;

    arb_init(lambda);
    arb_set_str(lambda, "1.1122", prec);
    arb_init(logQ);
    arb_init(r);
    arb_log_ui(logQ, 10000000000, prec);
    arb_div(r, lambda, logQ, prec);

    arb_t temp1;
    arb_init(temp1);
    arb_t temp2;
    arb_init(temp2);
    arb_t temp3;
    arb_init(temp3);
    arb_t pi;
    arb_init(pi);
    arb_const_pi(pi, prec);

    //Calculate C's
    arb_t C0;
    arb_init(C0);
    arb_set_str(C0, "2.97655", prec);
    arb_t logC[4];
    arb_init(logC[0]);
    arb_log(logC[0], C0, prec);
    for (int k=1; k<=3; k++) {
        // set P= 1/(4A)+5/(32A^2) = (32A^2-20A)/(128A^3) where A=(2^k)^2
        int A = (int)pow(2,2*k);
        arb_t P;
        arb_init(P);
        arb_t top;
        arb_init(top);
        arb_t bottom;
        arb_init(bottom);
        arb_set_ui(top, 32*A*A-20*A);
        arb_set_ui(bottom, 128*A*A*A);
        arb_div(P, top, bottom, prec);
        // set Pdiv = P/2^(k+2)
        arb_t Pdiv;
        arb_init(Pdiv);
        arb_div_ui(Pdiv, P, (int)pow(2,k+2), prec);

        // set Cterm = logC[k-1])/2
        arb_t Cterm;
        arb_init(Cterm);
        arb_div_ui(Cterm, logC[k-1], 2, prec);

        // set zetasum = sum_{j=0}^100 (-1)^j * \log\zeta(1+ 1/2^(k+1) + j/2^k)
        arb_t zetasum;
        arb_init(zetasum);
        arb_zero(zetasum);
        arb_t zetaterm;
        arb_init(zetaterm);
        arb_t zetainput;
        arb_init(zetainput);
        int parity = 1;
        for (int j=0; j<=1000; j++) {
            arb_set_ui(temp1, (int)pow(2,(k+1)));
            arb_inv(temp1, temp1, prec);
            arb_set_ui(temp2, j);
            arb_div_ui(temp2, temp2, (int)pow(2,k), prec);
            arb_add(zetainput, temp1, temp2, prec);
            arb_add_ui(zetainput, zetainput, 1, prec);
            arb_zeta(zetaterm, zetainput, prec);
            arb_log(zetaterm, zetaterm, prec);
            if (parity==-1) {
                arb_neg(zetaterm, zetaterm);
            }
            arb_add(zetasum, zetasum, zetaterm, prec);
            parity *= -1;
        }

        //set logC[k] = P/2^(k+2) + logC[k-1]/2 + zetasum
        arb_init(logC[k]);
        arb_add(logC[k], Pdiv, Cterm, prec);
        arb_add(logC[k], logC[k], zetasum, prec);

        printf("L(1-1/(2^(%d+1))) has constant ", k);
        arb_exp(temp1, logC[k], prec);
        printf(arb_get_str(temp1, 9,0));
        printf("\n");
    }

    //set fAB to be A/B
    arb_t f98;
    arb_init(f98);
    arb_set_ui(f98, 9);
    arb_div_ui(f98, f98, 8, prec);
    arb_t f54;
    arb_init(f54);
    arb_set_ui(f54, 5);
    arb_div_ui(f54, f54, 4, prec);
    arb_t f32;
    arb_init(f32);
    arb_set_ui(f32, 3);
    arb_div_ui(f32, f32, 2, prec);
    arb_t f12;
    arb_init(f12);
    arb_set_ui(f12, 1);
    arb_div_ui(f12, f12, 2, prec);
    arb_t f14;
    arb_init(f14);
    arb_set_ui(f14, 1);
    arb_div_ui(f14, f14, 4, prec);
    arb_t f18;
    arb_init(f18);
    arb_set_ui(f18, 1);
    arb_div_ui(f18, f18, 8, prec);
    arb_t f34;
    arb_init(f34);
    arb_set_ui(f34, 3);
    arb_div_ui(f34, f34, 4, prec);

    // set rAB to be r+A/B
    arb_t r18;
    arb_init(r18);
    arb_add(r18, r, f18, prec);
    arb_t r14;
    arb_init(r14);
    arb_add(r14, r, f14, prec);
    arb_t r12;
    arb_init(r12);
    arb_add(r12, r, f12, prec);
    arb_t r34;
    arb_init(r34);
    arb_add(r34, r, f34, prec);

    // sets r1 to r+7/8, the radius of the circle
    arb_t r1;
    arb_init(r1);
    arb_set_ui(temp1, 7);
    arb_div_ui(temp1, temp1, 8, prec);
    arb_add(r1, r, temp1, prec);

    // defines theta_i for i=1,2,3,4 as in the paper
    arb_t theta1;
    arb_init(theta1);
    arb_div(theta1, r18, r1, prec);
    arb_acos(theta1, theta1, prec);
    arb_t theta2;
    arb_init(theta2);
    arb_div(theta2, r14, r1, prec);
    arb_acos(theta2, theta2, prec);
    arb_t theta3;
    arb_init(theta3);
    arb_div(theta3, r12, r1, prec);
    arb_acos(theta3, theta3, prec);
    arb_t theta4;
    arb_init(theta4);
    arb_div(theta4, r34, r1, prec);
    arb_acos(theta4, theta4, prec);

    // sets sin_i = sin(theta_i), cos_i = cos(theta_i)
    arb_t sin1;
    arb_init(sin1);
    arb_t cos1;
    arb_init(cos1);
    arb_sin(sin1, theta1, prec);
    arb_cos(cos1, theta1, prec);

    arb_t sin2;
    arb_init(sin2);
    arb_t cos2;
    arb_init(cos2);
    arb_sin(sin2, theta2, prec);
    arb_cos(cos2, theta2, prec);

    arb_t sin3;
    arb_init(sin3);
    arb_t cos3;
    arb_init(cos3);
    arb_sin(sin3, theta3, prec);
    arb_cos(cos3, theta3, prec);

    arb_t sin4;
    arb_init(sin4);
    arb_t cos4;
    arb_init(cos4);
    arb_sin(sin4, theta4, prec);
    arb_cos(cos4, theta4, prec);

    // sets i_k and I_k to the integrals in the paper
    arb_t i1;
    arb_init(i1);
    arb_sub_ui(i1, sin1, 1, prec);
    arb_t i2;
    arb_init(i2);
    arb_sub(i2, sin2, sin1, prec);
    arb_t i3;
    arb_init(i3);
    arb_sub(i3, sin3, sin2, prec);
    arb_t i4;
    arb_init(i4);
    arb_sub(i4, sin4, sin3, prec);
    arb_t i5;
    arb_init(i5);
    arb_neg(i5, sin4);

    arb_t I1;
    arb_init(I1);
    arb_div_ui(temp1, pi, 2, prec);
    arb_sub(temp1, temp1, theta1, prec);
    arb_mul(temp2, sin1, cos1, prec);
    arb_sub(I1, temp1, temp2, prec);
    arb_div_ui(I1, I1, 2, prec);

    arb_t I2;
    arb_init(I2);
    arb_mul(temp1, sin1, cos1, prec);
    arb_mul(temp2, sin2, cos2, prec);
    arb_sub(temp2, temp1, temp2, prec);
    arb_sub(temp1, theta1, theta2, prec);
    arb_add(I2, temp1, temp2, prec);
    arb_div_ui(I2, I2, 2, prec);

    arb_t I3;
    arb_init(I3);
    arb_mul(temp1, sin2, cos2, prec);
    arb_mul(temp2, sin3, cos3, prec);
    arb_sub(temp2, temp1, temp2, prec);
    arb_sub(temp1, theta2, theta3, prec);
    arb_add(I3, temp1, temp2, prec);
    arb_div_ui(I3, I3, 2, prec);

    arb_t I4;
    arb_init(I4);
    arb_mul(temp1, sin3, cos3, prec);
    arb_mul(temp2, sin4, cos4, prec);
    arb_sub(temp2, temp1, temp2, prec);
    arb_sub(temp1, theta3, theta4, prec);
    arb_add(I4, temp1, temp2, prec);
    arb_div_ui(I4, I4, 2, prec);

    arb_t I5;
    arb_init(I5);
    arb_mul(temp1, sin4, cos4, prec);
    arb_add(I5, theta4, temp1, prec);
    arb_div_ui(I5, I5, 2, prec);

    //set phi to be the coefficient on logq
    arb_t phi;
    arb_init(phi);
    arb_div(phi, r, r1, prec);
    arb_mul(phi, phi, sin1, prec);
    arb_neg(phi, phi);
    arb_div(temp1, I1, r18, prec);
    arb_div_ui(temp1, temp1, 8, prec);
    arb_add(phi, phi, temp1, prec);
    arb_mul(temp2, sin1, cos1, prec);
    arb_add(temp2, temp2, theta1, prec);
    arb_div_ui(temp2, temp2, 2, prec);
    arb_add(phi, phi, temp2, prec);
    arb_div(phi, phi, pi, prec);

    //set R to be the coefficient on log(1+1/r)
    arb_t R;
    arb_init(R);
    arb_sub_ui(temp1, i1, 1, prec);
    arb_mul_ui(temp1, temp1, 2, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I1, r18, prec);
    arb_mul_ui(temp2, temp2, 2, prec);
    arb_div(temp2, temp2, pi, prec);
    arb_add(R, temp2, temp1, prec);
    arb_neg(R, R);

    // sets J0 equal to the J0 integral minus the term with \log(1+1/r)
    arb_t J0;
    arb_init(J0);
    arb_t sqrt_1; //equal to \sqrt{16r+15}
    arb_init(sqrt_1);
    arb_mul_ui(sqrt_1, r, 16, prec);
    arb_add_ui(sqrt_1, sqrt_1, 15, prec);
    arb_sqrt(sqrt_1, sqrt_1, prec);
    arb_t sqrt_2; //equal to \sqrt{16r+7}
    arb_init(sqrt_2);
    arb_mul_ui(sqrt_2, r, 16, prec);
    arb_add_ui(sqrt_2, sqrt_2, 7, prec);
    arb_sqrt(sqrt_2, sqrt_2, prec);
    arb_t sqrt7; //equal to \sqrt{7}
    arb_init(sqrt7);
    arb_sqrt_ui(sqrt7, 7, prec);

    arb_inv(temp1, sqrt_1, prec);
    arb_atan(temp1, temp1, prec);
    arb_mul(temp1, temp1, sqrt_1, prec);
    arb_mul_ui(temp1, temp1, 2, prec);

    arb_div(temp2, sqrt7, sqrt_2, prec);
    arb_atanh(temp2, temp2, prec);
    arb_mul(temp2, temp2, sqrt_2, prec);
    arb_mul(temp2, temp2, sqrt7, prec);
    arb_mul_ui(temp2, temp2, 2, prec);

    arb_mul_ui(J0, pi, 4, prec);
    arb_sub(J0, J0, temp1, prec);
    arb_sub(J0, J0, temp2, prec);
    arb_div(J0, J0, r1, prec);
    arb_div(J0, J0, r1, prec);
    arb_div_ui(J0, J0, 4, prec);
    arb_div(J0, J0, pi, prec);

    // calculates J1
    arb_t J1;
    arb_init(J1);
    arb_add(temp1, f98, r, prec);
    arb_mul(temp1, temp1, temp1, prec);
    arb_mul(temp2, r1, r1, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_log(temp1, temp1, prec);
    arb_div_ui(temp1, temp1, 16, prec);
    arb_mul_ui(temp2, logC[2], 2, prec);
    arb_add(J1, temp1, temp2, prec);
    arb_mul(J1, J1, I1, prec);
    arb_div(J1, J1, pi, prec);
    arb_div(J1, J1, r18, prec);

    // calculates J2
    arb_t J2;
    arb_init(J2);
    arb_t J2_1; // the first line in J2
    arb_init(J2_1);
    arb_t J2_2; // the second line in J2
    arb_init(J2_2);
    arb_add(temp1, f54, r, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta1, prec);
    arb_log(temp2, temp2, prec);
    arb_mul_ui(temp1, logC[1], 16, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r18, i2, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I2, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J2_1, temp1, temp3, prec);
    arb_add(temp1, f98, r, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta1, prec);
    arb_log(temp2, temp2, prec);
    arb_div_ui(temp2, temp2, 2, prec);
    arb_mul_ui(temp1, logC[2], 16, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r14, i2, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I2, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J2_2, temp1, temp3, prec);
    arb_neg(J2_2, J2_2);
    arb_add(J2, J2_1, J2_2, prec);

    // calculates J3
    arb_t J3;
    arb_init(J3);
    arb_t J3_1; // the first line in J3
    arb_init(J3_1);
    arb_t J3_2; // the second line in J3
    arb_init(J3_2);
    arb_add(temp1, f32, r, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta2, prec);
    arb_log(temp2, temp2, prec);
    arb_mul_ui(temp1, logC[0], 8, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r14, i3, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I3, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J3_1, temp1, temp3, prec);
    arb_add(temp1, f54, r, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta2, prec);
    arb_log(temp2, temp2, prec);
    arb_div_ui(temp2, temp2, 2, prec);
    arb_mul_ui(temp1, logC[1], 8, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r12, i3, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I3, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J3_2, temp1, temp3, prec);
    arb_neg(J3_2, J3_2);
    arb_add(J3, J3_1, J3_2, prec);

    // calculates J4
    arb_t J4;
    arb_init(J4);
    arb_t J4_1; // the first line in J4
    arb_init(J4_1);
    arb_t J4_2; // the second line in J4
    arb_init(J4_2);
    arb_t J4_3; // the third line in J4
    arb_init(J4_3);
    arb_add_ui(temp1, r, 2, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta3, prec);
    arb_log(temp2, temp2, prec);
    arb_mul_ui(temp1, pi, 2, prec);
    arb_log(temp1, temp1, prec);
    arb_mul_ui(temp1, temp1, 2, prec);
    arb_sub(temp3, temp2, temp1, prec);
    arb_mul(temp1, r12, i4, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I4, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J4_1, temp1, temp3, prec);
    arb_sub(temp1, f12, r, prec);
    norm(temp2, temp1, r1, theta4, prec);
    arb_log(temp2, temp2, prec);
    arb_mul_ui(temp1, logC[0], 8, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r34, i4, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I4, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J4_2, temp1, temp3, prec);
    arb_neg(J4_2, J4_2);
    arb_sub(temp1, f14, r, prec);
    norm(temp2, temp1, r1, theta4, prec);
    arb_log(temp2, temp2, prec);
    arb_div_ui(temp2, temp2, 2, prec);
    arb_mul_ui(temp1, logC[1], 8, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_mul(temp1, r12, i4, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I4, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J4_3, temp1, temp3, prec);
    arb_add(J4, J4_1, J4_2, prec);
    arb_add(J4, J4, J4_3, prec);


    // calculates J4
    arb_t J5;
    arb_init(J5);
    arb_t J5_1; // the first line in J5
    arb_init(J5_1);
    arb_t J5_2; // the second line in J5
    arb_init(J5_2);
    arb_t J5_3; // the third line in J5
    arb_init(J5_3);
    arb_add_ui(temp1, r, 2, prec);
    arb_neg(temp2, r1);
    norm(temp2, temp1, temp2, theta4, prec);
    arb_log(temp2, temp2, prec);
    arb_mul_ui(temp1, pi, 2, prec);
    arb_log(temp1, temp1, prec);
    arb_mul_ui(temp1, temp1, 2, prec);
    arb_sub(temp3, temp2, temp1, prec);
    arb_mul(temp1, r12, i5, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I5, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J5_1, temp1, temp3, prec);
    arb_log(temp2, f98, prec);
    arb_mul_ui(temp2, temp2, 2, prec);
    arb_mul_ui(temp1, logC[1], 16, prec);
    arb_add(temp3, temp1, temp2, prec);
    arb_div(temp1, i5, pi, prec);
    arb_div(temp2, I5, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J5_2, temp1, temp3, prec);
    arb_neg(J5_2, J5_2);
    arb_mul_ui(temp3, logC[2], 16, prec);
    arb_mul(temp1, r34, i5, prec);
    arb_div(temp1, temp1, pi, prec);
    arb_div(temp1, temp1, r1, prec);
    arb_div(temp2, I5, pi, prec);
    arb_add(temp1, temp1, temp2, prec);
    arb_mul(J5_3, temp1, temp3, prec);
    arb_add(J5, J5_1, J5_2, prec);
    arb_add(J5, J5, J5_3, prec);

    // calculates the O(1) terms from zeta
    arb_t zetaTerms;
    arb_init(zetaTerms);
    arb_t digamma;
    arb_init(digamma);
    arb_add_ui(temp1, r, 3, prec);
    arb_div_ui(temp1, temp1, 2, prec);
    arb_digamma(digamma, temp1, prec);
    arb_div_ui(digamma, digamma, 2, prec);
    arb_mul_ui(temp1, pi, 2, prec);
    arb_log(temp1, temp1, prec);
    arb_const_euler(temp2, prec);
    arb_div_ui(temp2, temp2, 2, prec);
    arb_sub(temp1, temp2, temp1, prec);
    arb_add_ui(temp1, temp1, 1, prec);
    arb_add(zetaTerms, temp1, digamma, prec);

    // add up all the O(1) contributions
    arb_t constant;
    arb_init(constant);
    arb_add(constant, J0, J1, prec);
    arb_add(constant, constant, J2, prec);
    arb_add(constant, constant, J3, prec);
    arb_add(constant, constant, J4, prec);
    arb_add(constant, constant, J5, prec);
    arb_add(constant, constant, zetaTerms, prec);

    // add up everything
    arb_t total; // total RHS
    arb_init(total);
    arb_t errorTerms; // RHS excluding the log(q) term
    arb_init(errorTerms);
    arb_mul(temp1, phi, logQ, prec);
    arb_inv(temp2, r, prec);
    arb_add_ui(temp2, temp2, 1, prec);
    arb_log(temp2, temp2, prec);
    arb_mul(temp2, temp2, R, prec);
    arb_add(errorTerms, temp2, constant, prec);
    arb_add(total, temp1, errorTerms, prec);

    // output values
    int decimalPlaces = 13;
    printf("phi: ");
    printf(arb_get_str(phi, decimalPlaces, 0));
    printf("\n");
    printf("R: ");
    printf(arb_get_str(R, decimalPlaces, 0));
    printf("\n");
    printf("total constant: ");
    printf(arb_get_str(constant, decimalPlaces, 0));
    printf("\n");
    printf("J0: ");
    printf(arb_get_str(J0, decimalPlaces, 0));
    printf("\n");
    printf("J1: ");
    printf(arb_get_str(J1, decimalPlaces, 0));
    printf("\n");
    printf("J2: ");
    printf(arb_get_str(J2, decimalPlaces, 0));
    printf("\n");
    printf("J3: ");
    printf(arb_get_str(J3, decimalPlaces, 0));
    printf("\n");
    printf("J4: ");
    printf(arb_get_str(J4, decimalPlaces, 0));
    printf("\n");
    printf("J5: ");
    printf(arb_get_str(J5, decimalPlaces, 0));
    printf("\n");
    printf("other terms: ");
    printf(arb_get_str(errorTerms, decimalPlaces, 0));
    printf("\n");
    printf("total: ");
    printf(arb_get_str(total, decimalPlaces, 0));
    printf("\n");


    return 0;
}

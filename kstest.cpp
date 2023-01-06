/**************************************************************************/
/*    Copyright (C) 2006 Romain Michalec                                  */
/*                                                                        */
/*    This library is free software; you can redistribute it and/or       */
/*    modify it under the terms of the GNU Lesser General Public          */
/*    License as published by the Free Software Foundation; either        */
/*    version 2.1 of the License, or (at your option) any later version.  */
/*                                                                        */
/*    This library is distributed in the hope that it will be useful,     */
/*    but WITHOUT ANY WARRANTY; without even the implied warranty of      */
/*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU   */
/*    Lesser General Public License for more details.                     */
/*                                                                        */
/*    You should have received a copy of the GNU Lesser General Public    */
/*    License along with this library; if not, write to the Free Software   */
/*    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,          */
/*    MA  02110-1301, USA                                                 */
/*                                                                        */
/**************************************************************************/

#include "kstest.h"
#include "KolmogorovSmirnovDist.hpp"
#include "Dsuite_utils.h"
#include <cmath>



/* KOLMOGOROV-SMIRNOV TEST OF HOMOGENEITY
 *
 * This long explaination provides you with everything you must know to
 * understand totally the Kolmogorov-Smirnov test of homogeneity.
 *
 * When x = (x_1, ..., x_n) are the observed values (the "realisation")
 * of a random sample X = (X_1, ..., X_n), we define the "sample cumulative
 * distribution function", or "empirical cumulative distribution function"
 * of X, and denote it X_emp_cdf(x), as the proportion of observed values
 * in X that are less or equal to x. In other words, if k observed values
 * in the sample are less than or equal to x, then X_emp_cdf(x) = k/n.
 *
 * Thus X_emp_cdf is a step function, defined on R and with values
 * from 0 to 1, with jumps of magnitude 1/n at each observation.
 *
 * The (non-empirical) cumulative distribution function of X, X_cdf,
 * defined on R, with values from 0 to 1, is X_cdf(x) = Pr (X_k <= x),
 * with k in 1:n (it doesn't matter which one since they all have the
 * same distribution).
 *
 * The empirical c.d.f is of course an approximation of the c.d.f:
 *
 *     X_emp_cdf(x) --> X_cdf(x)   when n --> infty
 *
 * In fact, a much stronger result, known as the Glivenko-Cantelli theorem,
 * states that when sample size is large, X_emp_cdf converges uniformly
 * to X_cdf:
 *
 *     D_n = sup |X_emp_cdf(x) - X_cdf(x)| --> 0   when n --> infty
 *
 * the sup function is over all real numbers (i.e. sup_{x \in R}).
 *
 * We consider a problem in which random samples are taken from two
 * populations: sample1, or X_1, ..., X_n, from population 1; sample2,
 * or Y_1, ..., Y_m, from population 2.
 * We must determine, on the basis of the observed values in the
 * samples, whether they come from the same distribution or not
 * (hence the name "test of homogeneity").
 *
 * The hypotheses to be tested are as follow:
 *
 *     H0: sample1 and sample2 are taken from distributions having
 *         the same distribution function
 *     H1: sample1 and sample2 are taken from distributions having
 *         different distribution functions
 *
 * To determine which hypothesis is the most likely, we must shape
 * a test statistic. Since the sample c.d.f.'s appeared to be estimators
 * of the non-empirical c.d.f.'s, and since c.d.f.'s are caracteristic of
 * the distribution of their random variable (i.e., if two random variables
 * have same c.d.f., then they also have same d.f.), it seems a good idea
 * to base a test on the difference between the two sample c.d.f.'s, and take
 * large sample sizes to approximate the difference between the corresponding
 * distribution functions.
 *
 * Hence we consider the statistic D_nm defined as follows:
 *
 *     D_nm = sup | X_emp_cdf(x) - Y_emp_cdf(x) |
 *
 * When the null hypothesis H0 is true, X_cdf and Y_cdf are identical functions
 * and we can easily deduce from the Glivenko-Cantelli theorem, by bounding
 * |X_emp_cdf(x) - Y_emp_cdf(x)|, that:
 *
 *     D_nm --> 0 when n,m --> infty
 *
 * It seems therefore reasonable to use a test of the form:
 *
 *     Reject H0 if D_nm > some critical value c
 *
 * but in fact, we will use the following test:
 *
 *     Reject H0 if sqrt(n*m/(n+m))*D_nm > c
 *
 * because of the following result, established in 1933 by Andreï Nikolaevitch
 * Kolmogorov and Vladimir Ivanovitch Smirnov:
 *
 *     For any t > 0,
 *     Pr (sqrt(n)*D_n <= t) -->
 *               1 - 2*Sum_{i=1:infty} (-1)^(i-1) exp(-2 i^2 t^2)
 *     when n --> infty
 *
 * a result that has been adapted to two-sample testing:
 *
 *     For any t > 0, if the null hypothesis H0 is true, then
 *     Pr (sqrt(n_approx)*D_nm <= t) -->
 *               1 - 2*Sum_{i=1:infty} (-1)^(i-1) exp(-2 i^2 t^2)
 *     when n,m --> infty and with n_approx = n*m/(n+m)
 *
 * We call "Kolmogorov's statistic" the statistic D_nm, and "limiting form
 * of Kolmogorov's statistic" the function
 *
 *     L(t) = 1 - 2*Sum_{i=1:infty} (-1)^(i-1) exp(-2 i^2 t^2)
 *
 * The result of Kolmogorov and Smirnov is of crucial importance, because
 * it enables the computation of Pr ( sqrt(nm/n+m)*D_nm <= t ) through the
 * (much easier) computation of the limiting form L(t), for large samples.
 *
 * As for Pearson chi-square test or Wilcoxon-Mann-Whitney ranks test,
 * we are interested in computing p-values for this test.
 *
 * For instance, let's suppose the observed values of sample1 and sample2
 * yield the value d for Kolmogorov's statistic D_nm, i.e. the value
 * sqrt(nm/n+m)*d for the statistic sqrt(nm/n+m)*D_nm. The very definition
 * of p as the significance level alpha corresponding to the critical value
 * c = sqrt(nm/n+m)*d reads:
 *
 *     p = Pr (H0 rejected | H0 true)
 *       = Pr (sqrt(nm/n+m)*D_nm > c = sqrt(nm/n+m)*d | H0 true)
 *       = Pr (D_nm > d | H0 true)
 *       = 1 - Pr (D_nm < d | H0 true)
 *
 * Moreover, we know that when H0 is true and sample sizes are large,
 *
 *     Pr (sqrt(nm/n+m)*D_nm < t) ~= L(t)
 *
 * i.e. with t = sqrt(nm/n+m)*d
 *
 *     Pr (D_nm < d) ~= L(sqrt(nm/n+m)*d)
 *
 * (we use ~= for "approximately equal to").
 *
 * Summary: once we have the value d of the test statistic D_nm (a value
 * that is really easy to obtain from the values in the samples),
 * we compute the p-value of the observed sample thanks to the formula
 * p = 1 - L(sqrt(nm/n+m)*d) (we stop the Sum in L when we think a
 * sufficient accuracy has been reached).
 *
 * In fact, in this implementation of Kolmogorov-Smirnov test of homogeneity,
 * we don't use the limiting form, but a small C procedure provided by
 * Marsaglia et al., K(int n,double d) (source code below), that computes
 * very quickly K(n,d) = Pr(D_n < d) with 13-15 digit accuracy -- much more
 * that what we need. The fact that it computes Pr(D_n < d) rather than
 * Pr(D_nm < d) is not a problem as we learned previously that both had the
 * same limiting forms, that samples sizes are large and that we do not
 * require a really high accuracy...
 *
 * Hence we have at last:
 *
 *     p = 1 - K(n_approx,d)    with  n_approx = nm/n+m
 *
 * Last thing to know: the differences between the one-sided and two-sided
 * versions of this test.
 *
 * "Reject H0 if sqrt(n*m/(n+m))*D_nm > c" is clearly one-sided; however,
 * because of the |.| in D_nm, the _meaning_ of the test is in fact two-sided.
 * Indeed, as D_nm = sup |X_emp_cdf(x) - Y_emp_cdf(x)|, rejecting H0 when
 * sqrt(nm/n+m)*D_nm > c means rejecting it when the "gap" between
 * X_emp_cdf(x) and Y_emp_cdf(x) is too large, in a direction or in the other,
 * i.e. X_emp_cdf(x) "above" or "under" Y_emp_cdf(x). Which translates
 * respectively into the X_k being globally smaller or larger than the Y_l.
 *
 * It is quite easy to make a two-sided test (with the meaning of a one-sided
 * test) from this one-sided test (with the meaning of a two-sided). This test
 * would be "Reject H0 if sqrt(nm/n+m)*D_nm != 0", or, better, the approximated
 * two-sided test: "Reject H_0 if |sqrt(nm/n+m)*D_nm| > epsilon", with epsilon
 * a very small number. However, this would not make much sense if D_nm was
 * still defined as sup |X_emp_cdf(x) - Y_emp_cdf(x)|; we should better
 * redefine it as D_nm = sup (X_emp_cdf(x) - Y_emp_cdf(x)).
 *
 * That would mean that, instead of rejecting H0 when the "gap" between
 * the c.d.f. X_emp_cdf(x) and Y_emp_cdf(x) is too large, in a direction or
 * in the other, we would reject it only if the "gap" is too large in the
 * sense of X_emp_cdf(x) being "above" Y_emp_cdf(x). Translation: when the
 * X_k are globally smaller than the Y_l. For instance, the X_k could be
 * the "old" traffic samples, and the Y_l the "new" ones, and we would
 * consider as anomalies only those that are more likely to be attacks
 * (increases in network traffic).
 *
 * As the p-value depends on the form of the test, we have to find another
 * formula for the p-value of this two-sided test, "Reject H_0 if
 * |sqrt(nm/n+m)*D_nm| > epsilon". Let's denote d the value of the
 * statistic D_nm yielded by the values in the sample, i.e. the value
 * of the statistic sqrt(nm/n+m)*D_nm is sqrt(nm/n+m)*d. Note that D_nm
 * is no longer Kolmogorov's statistic, but a modified Kolmogorov's statistic:
 * D_nm = sup (X_emp_cdf(x) - Y_emp_cdf(x)). As previously:
 *
 *     p = Pr (H0 rejected | H0 true)
 *       = Pr (|sqrt(nm/n+m)*D_nm| > sqrt(nm/n+m)*d | H0 true)
 *
 * Here we need a little drawing representing the distribution function of
 * the random variable sqrt(nm/n+m)*D_nm (let's assume it is something like
 * a regular curve from -infty to +infty, somewhat symetrical around 0 --
 * wishful thinking: there are no reasons it should be so):
 *
 *     p is the sum of the areas under the curve that are located further
 *     away from 0 than the value sqrt(nm/n+m)*d, i.e. from -infty to
 *     -sqrt(nm/n+m)*d and from sqrt(nm/n+m)*d to infty
 *
 * We suppose these areas are symetrical, each equals p/2 and we have:
 *
 *     1 - p/2 is the area under the curve from -infty to sqrt(nm/n+m)*d
 *
 * i.e.:
 *
 *     1 - p/2 = Pr (sqrt(nm/n+m)*D_nm < sqrt(nm/n+m)*d)
 *             = Pr (D_nm < d)
 *     p = 2 * (1 - Pr (D_nm < d))
 *
 * The difficulty resides in computing Pr (D_nm < d) now that D_nm is no
 * longer Kolmogorov's statistic... Neither the limiting form nor Marsaglia's
 * procedure are able to compute it, and although there exists a two-sample
 * version for Kolmogorov-Smirnov c.d.f. test, we are limited to the
 * one-sample version (unless we find one day a method to compute this
 * probability).
 */


/* The following three functions are copied from
 * G. Marsaglia, W. W. Tsang, J. Wang: "Evaluating  Kolmogorov's Distribution"
 *
 * The third one compute K(n,d) = Pr(D_n < d), where D_n is Kolmogorov's
 * goodness-of-fit measure for a sample c.d.f., n being the size of the
 * random sample: D_n = sup |X_emp_cdf(x) - X_cdf(x)| with the notations
 * of the previous text.
 *
 * The results correspond more or less to the tables found in
 * Hartung: "Statistik", 13rd Edition, pp. 521-523
 */
void mMultiply(double *A,double *B,double *C,int m)
{
    int i,j,k; double s;
    for(i=0;i<m;i++) for(j=0; j<m; j++)
    {s=0.; for(k=0;k<m;k++) s+=A[i*m+k]*B[k*m+j]; C[i*m+j]=s;}
}

void mPower(double *A,int eA,double *V,int *eV,int m,int n)
{
    double *B;int eB,i;
    if(n==1) {for(i=0;i<m*m;i++) V[i]=A[i];*eV=eA; return;}
    mPower(A,eA,V,eV,m,n/2);
    B=(double*)malloc((m*m)*sizeof(double));
    mMultiply(V,V,B,m); eB=2*(*eV);
    if(n%2==0){for(i=0;i<m*m;i++) V[i]=B[i]; *eV=eB;}
    else {mMultiply(A,B,V,m); *eV=eA+eB;}
    if(V[(m/2)*m+(m/2)]>1e140) {for(i=0;i<m*m;i++) V[i]=V[i]*1e-140;*eV+=140;}
    free(B);
}

double K(int n,double d)
{
   int k,m,i,j,g,eH,eQ;
   double h,s,*H,*Q;
    /* OMIT NEXT TWO LINES IF YOU REQUIRE >7 DIGIT ACCURACY IN THE RIGHT TAIL*/
s=d*d*n;
if(s>7.24||(s>3.76&&n>99) || n > 15000) return 1-2*exp(-(2.000071+.331/sqrt(n)+1.409/n)*s);
   k=(int)(n*d)+1;
   m=2*k-1;
   h=k-n*d;
   H=(double*)malloc((m*m)*sizeof(double));
   Q=(double*)malloc((m*m)*sizeof(double));
       for(i=0;i<m;i++)
         for(j=0;j<m;j++)
           if(i-j+1<0) H[i*m+j]=0;
          else     H[i*m+j]=1;
    for(i=0;i<m;i++)
    {
    H[i*m]-=pow(h,i+1);
    H[(m-1)*m+i]-=pow(h,(m-i));
    }
    H[(m-1)*m]+=(2*h-1>0?pow(2*h-1,m):0);
    for(i=0;i<m;i++)
    for(j=0;j<m;j++)
    if(i-j+1>0)
        for(g=1;g<=i-j+1;g++) H[i*m+j]/=g;
    eH=0;
    mPower(H,eH,Q,&eQ,m,n);
    s=Q[(k-1)*m+k-1];
    for(i=1;i<=n;i++)
    {
    s=s*i/n;
    if(s<1e-140){s*=1e140; eQ-=140;}
    }
    s*=pow(10.,eQ);
    std::cerr << "s: " << s << std::endl;
    
    free(H);
    free(Q);
    return s;
}


/* This function returns an approximated p-value of the Kolmogorov-Smirnov
 * c.d.f. test of homogeneity.
 *
 * The hypothesis to be tested are:
 *
 *     H0: the samples are drawn from populations having same distribution
 *     H1: they come from populations with different distributions
 *
 * The test statistic is (notations from the long explanation at the
 * beginning of this file):
 *
 *     D_nm = sup | X_emp_cdf(x) - Y_emp_cdf(x) |
 *
 * The p-value when the observed value for D_nm is d is computed through
 * (proof in the text at the beginning of this file):
 *
 *     p = 1 - Pr (D_nm < d)
 *       = 1 - K(n_approx,d)    with  n_approx = nm/n+m
 *
 * We reject H0 to significance level alpha if p-value < alpha,
 * i.e. considering the data we have, we reject H0 to any significance
 * level alpha > p-value.
 *
 * As explained at the end of the "opening text", Kolmogorov-Smirnov test
 * is a one-sided test (with a meaning of a two-sided), and the corresponding
 * two-sided test (with the meaning of a one-sided), although it exists and
 * makes sense, cannot be computed easily with our current means (Marsaglia's
 * procedure or the limiting form of D_nm).
 */
double ks_test (std::list<int64_t> sample1, std::list<int64_t> sample2,
        std::ostream& outfile, bool printDebug) {

  unsigned int n1, n2, n_approx;
    // sample sizes
  float d;
    // the value of Kolmogorov's statistic
    // for the particular values of sample1 and sample2
  int D, Dmin, Dmax, s;
    // used in computing this value d

  std::list<int64_t>::iterator it1, it2;

  // Determine sample sizes
  n1 = sample1.size();
  n2 = sample2.size();

  // Calculate a conservative n approximation
  n_approx = (unsigned) ceil(float(n1*n2)/(n1+n2));
  if (printDebug) outfile << "n_approx=" << n_approx << std::endl;

  // Sort samples
  sample1.sort(); //outfile << "sorted sample1: " << sample1 << std::endl;
  sample2.sort(); //outfile << "sorted sample2: " << sample2 << std::endl;

  // We divide the range 0..1 into n1*n2 intervals of equal size 1/(n1*n2).
  //
  // Each item in sample1 makes the sample c.d.f of sample1
  // jump by a step of n2 intervals.
  // Each item in sample2 makes the sample c.d.f of sample2
  // jump by a step of n1 intervals.
  //
  // For each item we compute D, related to the distance between the two
  // sample c.d.f., s_cdf_1 - s_cdf_2, by:
  //
  //    D/(n1*n2) = s_cdf_1 - s_cdf_2
  //
  // We want to determine:
  //
  //    Dmin/(n1*n2) = min [s_cdf_1 - s_cdf_2] <= 0
  //    Dmax/(n1*n2) = max [s_cdf_1 - s_cdf_2] >= 0
  //
  // And then the value of Kolmorogov's statistic D_n1n2 is just:
  //
  //    D_n1n2 = sup |s_cdf_1 - s_cdf_2|
  //           = max [ |Dmin/(n1*n2)| ; |Dmax/(n1*n2)| ]

  D = 0; Dmin = 0; Dmax = 0;
  it1 = sample1.begin();
  it2 = sample2.begin();

  while ( (it1 != sample1.end()) && (it2 != sample2.end()) ) {

    if (*it1 == *it2) {

        if (printDebug) outfile << *it1 << " tie!";
      // steps in both sample c.d.f., we need to perform all steps
      // in this point before comparing D to Dmin and Dmax

      s = *it1;
      // perform all steps in s_cdf_1 first
      do {
    D += n2;
    it1++;
      }
      while ( (*it1 == s) && (it1 != sample1.end()) );
      // perform all steps in s_cdf_2 now
      do {
    D -= n1;
    it2++;
      }
      while ( (*it2 == s) && (it2 != sample2.end()) );

      // now adapt Dmin, Dmax if necessary
      if (D > Dmax)
    Dmax = D;
      else if (D < Dmin)
    Dmin = D;

    }

    else if (*it1 < *it2) {

    if (printDebug) outfile << *it1;
      // step in s_cdf_1, increase D by n2
      D += n2;
      it1++;

      if (D > Dmax)
    Dmax = D;

    }

    else {

      if (printDebug) outfile << *it2;
      // step in F2, decrease D by n1
      D -= n1;
      it2++;

      if (D < Dmin)
    Dmin = D;

    }

      if (printDebug) outfile << " D=" << D << " Dmin=" << Dmin << " Dmax=" << Dmax << std::endl;

  }

  // For two-sided test, take D = max (|Dmax|, |Dmin|) and compute
  // the value d of Kolmogorov's statistic (two-sided only)

  if (-Dmin > Dmax)
    D = -Dmin;
  else
    D = Dmax;

  // Hence the observed value of Kolmogorov's statistic:
  d = float(D)/(n1*n2);

  // Return p-value
  return 1 - K(n_approx,d);

}


/*
 D+ = max{(i/N)-Ri}, 1<=i<=N
 D- = max{(Ri-((i-1)/N)}, 1<=i<=N

 Step4: Compute calculated D:
 D= max(D+, D-):
 
 
 
 */

double ks_test_of_uniformity(std::vector<double> sampleVect0to1, std::ostream& outfile, bool printDebug) {

    unsigned int N = sampleVect0to1.size();
    
    double d, Dplusmax, Dminusmax;
    std::vector<double> DplusVals(N,0.0);
    std::vector<double> DminusVals(N,0.0);
    
    for (int i = 0; i < N; i++) {
        int j = i+1;
        double ratio = (double)j/N;
        double ratiominus = (double)i/N;
        DplusVals[i] = ratio - sampleVect0to1[i];
        DminusVals[i] = sampleVect0to1[i] - ratiominus;
    }
    
    Dplusmax = *max_element(DplusVals.begin(), DplusVals.end());
    Dminusmax = *max_element(DminusVals.begin(), DminusVals.end());
    
    
  //  std::cerr << "Dplusmax: " << Dplusmax << std::endl;
  //  std::cerr << "Dminusmax: " << Dminusmax << std::endl;
    
  //  print_vector(DplusVals, std::cerr, ',');
  //  print_vector(DminusVals, std::cerr, ',');
    
    
    if (Dplusmax > Dminusmax) {
        d = Dplusmax;
    } else {
        d = Dminusmax;
    }
    
   // std::cerr << "d: " << d << std::endl;
    
    // Return p-value
    return 1 - KScdf(N,d);

}

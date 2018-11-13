/*
Copyright 2018 Oleg Skoromnik

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

#include <config.h>
#include <effz_lib/effz_second_order_generating_int.h>
#include <effz_lib/effz_spec_func.h>

#include <cmath>
#include <iostream>

namespace effz{
	namespace second_order{
		namespace single{
			inline double der_NLpZpowMP(double z,
					unsigned int n, double lambda,
					int p, unsigned int i)
			{
				return binpow(n,i)
					* binpow((static_cast<double>(n) * lambda + z),p+i)
					* factorial_power(-p,i);
			}

			double der_LpL1powQ_div_NLpZpowP(double z,
					unsigned int n, double lambda, double lambda1,
					int q, int p, unsigned int i)
			{
				auto term = [&](unsigned int r){
					return binomial(i,r)
						* binpow(lambda+lambda1, q-(i-r))
						* factorial_power(q, i-r)
						* der_NLpZpowMP(z,n,lambda,p,r);
				};
				double res = 0.;
				for(unsigned int r = 0; r <= i; ++r){
					res += term(r);
				}
				return res;
			}

			double der_NLmZpowQ_div_NLpZpowP(double z,
					unsigned int n, double lambda,
					int q, int p, unsigned int i)
			{
				auto term = [&](unsigned int s){
					return binomial(i,s)
						* factorial_power(q,i-s)
						* factorial_power(-p,s)
						* binpow(static_cast<double>(n)*lambda-z,q-(i-s))
						/ binpow(static_cast<double>(n)*lambda+z,p+s);
				};
				double res = 0.;
				for(unsigned int s = 0; s <= i; ++s){
					res += term(s);
				}
				return binpow(n,i)*res;
			}

			double generating_general_term_a(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return -static_cast<double>((n+l+1))
					/ static_cast<double>(2*n)
					* binpow(4.*z*z,2*l+2)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l,n+l+2,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l,n+l+2,j);
			}

			double generating_general_term_b(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return (poly_gamma(n+l+1) - poly_gamma(n-l)
						- static_cast<double>(2*l+3)
						/ static_cast<double>(2*n))
					* binpow(4.*z*z,2*l+2)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l-1,n+l+1,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-1,n+l+1,j);
			}

			double generating_general_term_c(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				if(n-l-2 < 0) return 0.;

				auto der1 = [&](int k){
					return der_NLmZpowQ_div_NLpZpowP(z,n,lambda,
							k,k+2*l+2,i)
						* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,
								k,k+2*l+2,j);
				};
				double res = 0.;
				for(int k = 0; k <= n-l-2; ++k){
					res += (der1(k) / static_cast<double>(n-l-1-k));
				}
				return binpow(4.*z*z, 2*l+2) * res;
			}

			double generating_general_term_d(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return static_cast<double>(n+l+1)
					/ static_cast<double>(2*n)
					* binpow(4.*z*z,2*l+3)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l-1,n+l+2,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-1,n+l+2,j);
			}

			double generating_general_term_e(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return static_cast<double>(n-l-1)
					/ static_cast<double>(2*n)
					* binpow(4.*z*z,2*l+2)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l-2,n+l,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-2,n+l,j);
			}

			double generating_general_term_f(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return -static_cast<double>(n-l-1)
					/ static_cast<double>(2*n)
					* binpow(4.*z*z,2*l+3)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l-2,n+l+1,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-2,n+l+1,j);
			}

			double generating_general_term_g(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				double nLmZ = static_cast<double>(n)*lambda - z;
				double nLpZ = static_cast<double>(n)*lambda + z;
				double nL1mZ = static_cast<double>(n)*lambda1 - z;
				double nL1pZ = static_cast<double>(n)*lambda1 + z;
				double LpL1 = lambda+lambda1;
				double fourZZ_2lp2 = binpow(4.*z*z,2*l+2);

				if(i==0 && j==0){
					return fourZZ_2lp2 * binpow(nLmZ*nL1mZ, n-l-1)
						/ binpow(nLpZ*nL1pZ, n+l+1)
						* std::log(2.*static_cast<double>(n)
								* z * (lambda+lambda1)
								/ (nLpZ * nL1pZ));
				}

				auto der1_logNLpZfrac = [&n,&z](double L, int q, int p,
						unsigned int ip){
					auto term = [&](int s){
						if(s < 0) return 0.;
						int minus1_pow = (ip-s+1)&1 ? -1 : 1;
						return binomial(ip,s) * minus1_pow
							* binpow(n/(static_cast<double>(n)*L+z),ip-s)
							* factorial(ip-s-1)
							* der_NLmZpowQ_div_NLpZpowP(z,n,L,q,p,s);
					};
					double a = std::log(static_cast<double>(n)*L+z)
						* der_NLmZpowQ_div_NLpZpowP(z,n,L,q,p,ip);
					if(ip-1 < 0){
						return a;
					}
					for(int s = 0; s <= ip-1; ++s){
						a += term(s);
					}
					return a;
				};

				auto der1_log2NZLpL1_and_logNLpZ_frac
					= [&n,&z](double L, double L1,
						int q, int p, unsigned int ip){
					auto term = [&](int s){
						if(s < 0) return std::make_tuple(0.,0.);
						int minus1_pow = (ip-s+1)&1 ? -1 : 1;
						double factor = binomial(ip,s) * minus1_pow
							* factorial(ip-s-1)
							* der_NLmZpowQ_div_NLpZpowP(z,n,L,q,p,s);
						double der_log2NZ = binpow(L+L1, s-ip);
						double der_logNLpZ
							= binpow(n/(static_cast<double>(n)*L+z),ip-s);

						return std::make_tuple(
								factor * der_log2NZ,
								factor * der_logNLpZ
								);
					};
					double d = der_NLmZpowQ_div_NLpZpowP(z,n,L,q,p,ip);
					double log2NZLpL1_d = std::log(2.
							* static_cast<double>(n) * z * (L+L1)) * d;
					double logNLpZ_d
						= std::log(static_cast<double>(n) * L + z) * d;

					if(ip-1 < 0){
						return std::make_tuple(log2NZLpL1_d, logNLpZ_d);
					}
					for(int s = 0; s <= ip-1; ++s){
						auto tuple = term(s);
						log2NZLpL1_d += std::get<0>(tuple);
						logNLpZ_d += std::get<1>(tuple);
					}
					return std::make_tuple(log2NZLpL1_d, logNLpZ_d);
				};
				auto der1 = [&z,&n,&l,&der1_log2NZLpL1_and_logNLpZ_frac]
					(double L, double L1,
							 unsigned int ip){
					double NL1mZ = static_cast<double>(n)*L1 - z;
					double NL1pZ = static_cast<double>(n)*L1 + z;
					double prefactor = binpow(NL1mZ, n-l-1)
						/ binpow(NL1pZ, n+l+1);
					auto tuple = der1_log2NZLpL1_and_logNLpZ_frac(L,L1,
							n-l-1,n+l+1,ip);
					double der_log2NZLpL1_frac = std::get<0>(tuple);
					double der_logNLpZ_frac = std::get<1>(tuple);
					return prefactor * (der_log2NZLpL1_frac
							- der_logNLpZ_frac
							- std::log(NL1pZ)*der_NLmZpowQ_div_NLpZpowP(z,
								n,L,n-l-1,n+l+1,ip)
							);
				};

				if(i==0 && j!= 0)
					return fourZZ_2lp2 * der1(lambda1,lambda,j);
				if(j==0 && i!= 0)
					return fourZZ_2lp2 * der1(lambda,lambda1,i);

				auto term = [&](int s, int r){
					double minus1_pow = (1+i-s)&1 ? -1 : 1;
					return binomial(i,s) * minus1_pow
						* factorial(i-s-1)
						* der_NLmZpowQ_div_NLpZpowP(z,n,lambda,n-l-1,
								n+l+1,s)
						* binomial(j,r) * binpow(LpL1,-i-j+s+r)
						* factorial_power(s-i,j-r)
						* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-1,
								n+l+1,r);
				};
				double sum = 0.;
				if(i-1 >= 0){
					for(int s = 0; s <= i-1; ++s){
						for(int r = 0; r <= j; ++r){
							sum += term(s,r);
						}
					}
				}
				auto tuple = der1_log2NZLpL1_and_logNLpZ_frac(lambda1,
					  lambda, n-l-1, n+l+1, j); 
				double der_log2NZLpL1_frac_j = std::get<0>(tuple);
				double der_logNLpZ_frac_j = std::get<1>(tuple);
				auto der_same = der_NLmZpowQ_div_NLpZpowP(z,n,lambda,
						n-l-1,n+l+1,i);

				sum += der_log2NZLpL1_frac_j * der_same;
				sum -= der1_logNLpZfrac(lambda,n-l-1,n+l+1,i)
					* der_NLmZpowQ_div_NLpZpowP(z,n,lambda1,n-l-1,n+l+1,j);
				sum -= der_logNLpZ_frac_j * der_same;

				return fourZZ_2lp2 * sum;
			}

			double der_generating_general_z(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				double n_div_2z = static_cast<double>(n) / (2.*z);
				double prefactor = -2. * binpow(n_div_2z, 2*l+3)
					* factorial(n+l) / factorial(n-l-1);
				std::cout << prefactor
					* generating_general_term_a(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_b(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_c(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_d(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_e(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_f(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_g(z,n,l,lambda,lambda1,i,0)
					<< "\n";
				std::cout << prefactor
					* generating_general_term_g(z,n,l,lambda,lambda1,i,j)
					<< "\n";
				return prefactor * (
					generating_general_term_a(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_b(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_c(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_d(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_e(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_f(z,n,l,lambda,lambda1,i,j)
					+ generating_general_term_g(z,n,l,lambda,lambda1,i,j)
						);
			}
		} /* end namespace single */
	} /* end namespace second_order */
} /* end namespace effz */

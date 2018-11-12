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

#ifndef EFFZ_SPEC_FUNC_H
#define EFFZ_SPEC_FUNC_H

#include <complex>
#include <array>
#include <limits>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_psi.h>


namespace effz {
	/*
	 *abbreviations
	 *    h_*** --- hydrogen_***, i.e. charge = 1
	 *    h_l_*** --- hydrogen_like_***, i.e. charge != 1
	 *    ***_gsl function with a signature of a gsl one
	 */
	inline double factorial_power(int x, int n)
	{
		int res = 1;
		for(int i = 0; i < n; ++i){
			res *= (x-i);
		}
		return static_cast<double>(res);
	}

	template<typename T>
	inline T binpow(T a, int n)
	{
		T res = static_cast<T>(1);
		if(a == static_cast<T>(0) && n == 0) return res;

		int pow = n;
		T base = a;
		if(n < 0){
			if(a == 0) return std::numeric_limits<T>::infinity();
			pow = -n;
			base = 1./a;
		}

		while(pow){
			if(pow&1)
				res *= base;
			base *= base;
			pow >>= 1;
		}
		return res;
	}

	inline double binomial(unsigned int n, unsigned int m)
	{
		return gsl_sf_choose(n,m);
	}

	inline double poly_gamma(unsigned n)
	{
		return gsl_sf_psi_int(n);
	}

	inline double factorial(unsigned int n)
	{
		return gsl_sf_fact(n);
	}

	double h_rnl(const int n, const int l, const double r);

	struct h_rnl_params{
		int n;
		int l;
	};
	double h_rnl_gsl(double r, void *h_rnl_params_);

	double h_l_rnl(
			const double z,
			const int n,
			const int l,
			const double r
			);

	struct h_l_rnl_params{
		double z;
		int n;
		int l;
	};
	double h_l_rnl_gsl(double r, void *h_l_rnl_params_);

#ifndef EFFZ_ZEROTH_ORDER
	std::complex<double> green_coulomb_h_rad(
			const std::complex<double> E,
			const int l,
			const std::array<double, 2> &r
			);

	struct gf_h_params{
		std::complex<double> E;
		int l;
	};
	std::complex<double> green_coulomb_h_rad_gsl(
			double *r,
			size_t dim,
			void *gf_h_params_
			);

	std::complex<double> green_coulomb_rad(
			const double z,
			const std::complex<double> E,
			const int l,
			const std::array<double, 2> &r
			);

	struct gf_h_l_params{
		double z;
		std::complex<double> E;
		int l;
	};
	std::complex<double> green_coulomb_rad_gsl(
			double *r,
			size_t dim,
			void *gf_h_l_params_
			);
#endif

	inline double h_E(int n){
		double nn = (double)n;
		return - 1. / (2. * nn * nn);
	}

	inline double h_l_E(double z, int n){
		double nn = (double)n;
		return - z * z / (2. * nn * nn);
	}

	double three_j_symbol(const std::array<std::array<int,3>,2> &jm);

	std::complex<double> sph_harm_y(
			const int l,
			const int m,
			const double theta,
			const double phi
			);


} /*end namespace effz*/

#endif /* EFFZ_SPEC_FUNC_H */

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

#ifndef EFFZ_ZEROTH_ORDER_H
#define EFFZ_ZEROTH_ORDER_H

#include "effz_typedefs.h"

#ifdef __cplusplus

namespace effz{
	namespace zeroth_order{
		/*
		 *In this file "par" means parallel function
		 */
		double three_j_prod_direct(
				const int l,
				const int m,
				const int l1,
				const int m1,
				const int k);

		double i_direct(
				const int n,
				const int l,
				const int n1,
				const int l1,
				const int k);

		double three_j_prod_exchange(
				const int l,
				const int m,
				const int l1,
				const int m1,
				const int k);

		double i_exchange(
				const int n,
				const int l,
				const int n1,
				const int l1,
				const int k);

		double v_direct(
				const int n,
				const int l,
				const int m,
				const int n1,
				const int l1,
				const int m1);

		double v_exchange(
				const int n,
				const int l,
				const int m,
				const int n1,
				const int l1,
				const int m1);

		double v_direct_total(const occ_nums_array &g);

		double v_exchange_total(const occ_nums_array &g);

		double v_direct_total_par(const occ_nums_array &g);

		double v_exchange_total_par(const occ_nums_array &g);

		double v_total(const occ_nums_array &g);

		double v_total_par(const occ_nums_array &g);

		double a(const occ_nums_array &g);

		double z_star_0th(double z,
				const occ_nums_array &g);

		double e_0th(double z,const occ_nums_array &g);

		double z_star_0th_par(double z, const occ_nums_array &g);

		double e_0th_par(double z,const occ_nums_array &g);

		std::tuple<double,double> z_star_and_e_0th_par(double z,
				const occ_nums_array &g);

		class density_0th
		{
			private:
				double z;
				occ_nums_array occ_nums;
			public:
				density_0th(const double z,
						const occ_nums_array &occ_nums);
				double operator()(
						const double r,
						const double theta,
						const double phi) const;
		};
	} /*end namespace zeroth_order*/

} /*end namespace effz*/
#endif
#ifdef __cplusplus
extern "C" {
#endif
	double effz_three_j_prod_direct(
			const int l,
			const int m,
			const int l1,
			const int m1,
			const int k);

	double effz_i_direct(
			const int n,
			const int l,
			const int n1,
			const int l1,
			const int k);

	double effz_three_j_prod_exchange(
			const int l,
			const int m,
			const int l1,
			const int m1,
			const int k);

	double effz_i_exchange(
			const int n,
			const int l,
			const int n1,
			const int l1,
			const int k);

	double effz_v_direct(
			const int n,
			const int l,
			const int m,
			const int n1,
			const int l1,
			const int m1);

	double effz_v_exchange(
			const int n,
			const int l,
			const int m,
			const int n1,
			const int l1,
			const int m1);

	double effz_v_direct_total(const effz_occ_num_t *g, size_t dim);

	double effz_v_exchange_total(const effz_occ_num_t *g, size_t dim);

	double effz_v_direct_total_par(const effz_occ_num_t *g, size_t dim);

	double effz_v_exchange_total_par(const effz_occ_num_t *g, size_t dim);

	double effz_v_total(const effz_occ_num_t *g, size_t dim);

	double effz_v_total_par(const effz_occ_num_t *g, size_t dim);

	double effz_a(const effz_occ_num_t *g, size_t dim);

	double effz_z_star_0th(double z, const effz_occ_num_t *g, size_t dim);

	double effz_e_0th(double z, const effz_occ_num_t *g, size_t dim);

	double effz_z_star_0th_par(double z,
			const effz_occ_num_t *g, size_t dim);

	double effz_e_0th_par(double z, const effz_occ_num_t *g, size_t dim);

	/*
	 *density_0th class start
	 */
	typedef void* effz_density_0th_t;
	effz_density_0th_t effz_density_0th_new(double z,
			const effz_occ_num_t *nums, size_t dim);

	void effz_density_0th_delete(effz_density_0th_t rho);

	double effz_density_0th_at(const effz_density_0th_t rho,
			double r, double theta, double phi);
	/*
	 *density_0th class end
	 */


#ifdef __cplusplus
}
#endif
#endif /* EFFZ_ZEROTH_ORDER_H */

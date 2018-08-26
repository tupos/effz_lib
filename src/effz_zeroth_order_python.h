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

#ifndef EFFZ_ZEROTH_ORDER_PYTHON_H
#define EFFZ_ZEROTH_ORDER_PYTHON_H

#include <Python.h>

#include <effz_lib/effz_typedefs.h>

#ifdef __cplusplus
#include <vector>
#include <array>
#include <string>
#include <memory>
#include <functional>


namespace effz{
	namespace zeroth_order{

		PyObject* computed_rho_h_l(const occ_nums_array &g);
		PyObject* computed_rho_h_l_fourier(const occ_nums_array &g);
		PyObject* computed_asf_h_l(const occ_nums_array &g);
		void print_rho_h_l(const occ_nums_array &g);
		void print_rho_h_l_fourier(const occ_nums_array &g);
		void print_asf_h_l(const occ_nums_array &g);

		class symbolic_density{
			public:
				symbolic_density(const occ_nums_array &g);
				std::string get_density_latex_str();
				std::string get_density_pretty_str();
			private:
				const occ_nums_array g;
				std::unique_ptr<PyObject,std::function<void(PyObject*)>>
					density_ptr;
				std::wstring density_latex_str;
				std::wstring density_pretty_str;
		};

		class symbolic_asf{
			public:
				symbolic_asf(const occ_nums_array &g);
				std::string get_asf_latex_str();
				std::string get_asf_pretty_str();
			private:
				const occ_nums_array g;
				std::unique_ptr<PyObject,std::function<void(PyObject*)>>
					asf_ptr;
				std::wstring asf_latex_str;
				std::wstring asf_pretty_str;
		};

	} /* end namespace zeroth_order */
} /* end namespace effz */
#endif
#ifdef __cplusplus
extern "C" {
#endif
	PyObject* effz_computed_rho_h_l(const effz_occ_num_t *g, size_t dim);

	PyObject* effz_computed_rho_h_l_fourier(const effz_occ_num_t *g,
			size_t dim);

	PyObject* effz_computed_asf_h_l(const effz_occ_num_t *g, size_t dim);

	void effz_print_rho_h_l(const effz_occ_num_t *g, size_t dim);

	void effz_print_rho_h_l_fourier(const effz_occ_num_t *g, size_t dim);

	void effz_print_asf_h_l(const effz_occ_num_t *g, size_t dim);

	/*
	 *symbolic density class start
	 */
	typedef void* effz_symbolic_density_0th_t;
	effz_symbolic_density_0th_t effz_symbolic_density_0th_new(
			const effz_occ_num_t *g, size_t dim);

	void effz_symbolic_density_0th_delete(effz_symbolic_density_0th_t rho);

	int effz_symbolic_density_0th_get_pretty_str(
			const effz_symbolic_density_0th_t rho, char** str, size_t *len);

	int effz_symbolic_density_0th_get_latex_str(
			const effz_symbolic_density_0th_t rho, char** str, size_t *len);
	/*
	 *symbolic density class end
	 */

	/*
	 *symbolic asf class start
	 */
	typedef void* effz_symbolic_asf_0th_t;
	effz_symbolic_asf_0th_t effz_symbolic_asf_0th_new(
			const effz_occ_num_t *g, size_t dim);

	void effz_symbolic_asf_0th_delete(effz_symbolic_asf_0th_t rho);

	int effz_symbolic_asf_0th_get_pretty_str(
			const effz_symbolic_asf_0th_t asf, char** str, size_t *len);

	int effz_symbolic_asf_0th_get_latex_str(
			const effz_symbolic_asf_0th_t asf, char** str, size_t *len);
	/*
	 *symbolic asf class end
	 */
#ifdef __cplusplus
}
#endif
#endif /* EFFZ_ZEROTH_ORDER_PYTHON_H */

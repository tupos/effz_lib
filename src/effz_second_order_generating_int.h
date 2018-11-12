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

#ifndef EFFZ_SECOND_ORDER_GENERATING_INT_H
#define EFFZ_SECOND_ORDER_GENERATING_INT_H



namespace effz{
	namespace second_order{
		namespace single{
			double der_generating_general_z(double z, unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j);

			inline double der_generating_general(unsigned int n,
					unsigned int l, double lambda, double lambda1,
					unsigned int i, unsigned int j)
			{
				return der_generating_general_z(1,n,l,lambda,
						lambda1,i,j);
			}

		} /* end namespace single */
	} /* end namespace second_order */
} /* end namespace effz */

#endif /* EFFZ_SECOND_ORDER_GENERATING_INT_H */

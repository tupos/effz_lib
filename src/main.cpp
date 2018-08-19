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

#include <Python.h>
#include <config.h>

#include "effz_config.h"

#include <cstdlib>
#include <csignal>
#include <iostream>

using namespace eff_z;


int main(int argc, char *argv[]) try {

	config::shared_config().check_dirs();

	putenv((char*)"PYTHONDONTWRITEBYTECODE=1");
	Py_Initialize();
	PyRun_SimpleString("import sys\nimport os");
	PyRun_SimpleString(config::shared_config()
			.get_python_path_cmd().c_str());

#ifdef DEBUG
	PyRun_SimpleString("print(sys.path)");
#endif

	std::cout << "Hello world\n";

#ifdef EFFZ_PYTHON_VERSION_SIX
	if(Py_FinalizeEx() < 0){
		return 120;
	}
#elif defined(EFFZ_PYTHON_VERSION_FIVE)\
	|| defined(EFFZ_PYTHON_VERSION_FOUR)
	Py_Finalize();
#endif

	return EXIT_SUCCESS;

} catch(const std::exception &e){
	std::cerr << e.what() << "\n";
	return EXIT_FAILURE;
} catch(...){
	return EXIT_FAILURE;
}

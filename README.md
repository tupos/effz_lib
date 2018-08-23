## effz_lib - Library for the effective charge program

This is effz_lib, a library for the acompaniying effz program.

Contents
--------

- [Info](#info)
- [Contact](#contact)
- [License](#license)

Info
-----

The library provides `c` and `c++` api for the computation of the observables characteristics. 

The dependencies are described in the [effz][effz] program on github.

If the library is used in scientific publication we ask to provide the corresponding citation for the article 

	O D Skoromnik et al 2017 J. Phys. B: At. Mol. Opt. Phys. 50 245007

During install of the library there is a usefull option. This is to specify the `--with-effz-home` variable in `configure`. This is a directory in which the python code and databases are stored. By default they are stored under `$HOME/.config/effz_lib`. So the option `--with-effz-home=prefix` overrides the standard path `$HOME/.config` to the one provided by `prefix`.


Contact
-------

If you have questions about the program or need help, please write me
an email <olegskor@gmail.com>.

If you have bug reports or feature suggestions, please use the [issue
tracker][issueTracker].

The latest version of the program is available at
<https://github.com/tupos/effz_lib>.

The author's homepage is <https://github.com/tupos>.

License
-------

effz_lib is free software, you can redistribute it and/or modify it under
the terms of [Apache License 2.0][apache].

Unless required by applicable law or agreed to in writing, Licensor
provides the Work (and each Contributor provides its Contributions) on
an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either
express or implied, including, without limitation, any warranties or
conditions of TITLE, NON-INFRINGEMENT, MERCHANTABILITY, or FITNESS FOR
A PARTICULAR PURPOSE.

© 2018 effz authors

[apache]: https://www.apache.org/licenses/LICENSE-2.0
[arxiv]: https://arxiv.org/abs/1701.04800
[jphysb]: http://iopscience.iop.org/article/10.1088/1361-6455/aa92e6/meta
[issueTracker]: https://github.com/tupos/effz/issues
[effz]: https://github.com/tupos/effz

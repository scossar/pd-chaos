# PD Chaos

The beginnings of a library of Pure Data externals for various chaotic systems:

- Chua's circuit
- double pendulum
- Duffing oscillator
- Henon map
- Lorenz Attractor
- Logistic map
- Rossler attractor
- Tent map

The objects are usable, but should be considered to be under development for now - at this point the documentation is reading the source code.

## Building the externals

I've been using [pd-lib-builder](https://github.com/pure-data/pd-lib-builder) helper makefile to build the externals. To use it, clone the pd-chaos repo, then pull in the pd-lib-builder code with:

```bash
git submodule update --init --recursive
```
You should then be able to build the code with:

```bash
make install pdincludepath=/<path_to>/pure-data/src/ objectsdir=/<path_to>/pd/externals
```


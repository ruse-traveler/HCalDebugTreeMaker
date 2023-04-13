This repository contains a simple JANA plugin to prepare a tree of Barrel Hadronic Calorimeter cells (i.e. tiles) which can be used to train a ML clusterizer, `JBarrelHCalTreeMaker`. (Derived from code by Frederike Bock. Thanks!!)

### JBarrelHCalTreeMaker Usage
```
# after compiling EICrecon, do:
eicmkplugin.py JBarrelHCalTreeMaker
cp JBarrelHCalTreeMakerProcessor.* $EICrecon_ROOT/JBarrelHCalTreeMaker/
cmake -S JBarrelHCalTreeMaker -B JBarrelHCalTreeMaker/build
cmake --build JBarrelHCalTreeMaker/build --target install
eicrecon -Pplugins=JBarrelHCalTreeMaker <input edm4hep file>
```

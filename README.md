# switch_cell_behaviour

Switch behaviour of a cell based on its type. Here a single biomodule (`GeneralModule`) is used, that runs the simulation 
step based on the `cell_type`, either `type1behaviour()` or `type2behaviour()` is called in the modules' `Run()`routine.

Simulation starts with `cell_type=1` and is switched to `cell_type=2` when the  intracellular substance quantity `counter_sub_quant` falls below a predefind threshold.


##### Issues

Yet compilation fails with error:
```
/switch_cell_behaviour/src/switch_cell_behaviour.h:37:66: error: no matching constructor for initialization of 'bdm::Cell'
        explicit MyCell(const std::array<double, 3>& position) : Base(position) {
```
-----------------
##### BioDynaMo
```
$ biodynamo --version
v0.1.0-374-g59d6ef77
```

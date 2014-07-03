spline_ops
==========

Submodule with various operations on some common splines.

**STATUS**
==========

More or less functional, though the comments and docs could stand some improvement.

A few basic splines are unimplemented (TCB, uniform B splines), but not really high-priority.

It could stand a helper module, I imagine, since putting all the pieces together (say, integrands and arc lengths)
is rather non-obvious.

**DEPENDENCIES**
================

Tested on Lua 5.1.

The `arc_length` module depends on some [integrators](https://github.com/ggcrunchy/corona-sdk-snippets/blob/master/number_ops/integrators.lua)
support, which is in the rather unstable `number_ops` directory (logically it's four or five submodules, I think). It is very likely this will
move, say to a "solvers" submodule, but I'm still thinking about it.
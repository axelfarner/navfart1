# flypath3d_v3

This directory is a backward-compatible update of the original FlyPath3D
package by Witold Bużantowicz.

## Compatibility update

MATLAB R2025a introduced a new graphics renderer. With dense triangular
models, patch edges can be joined incorrectly and appear as long radiating
lines. Version 3 draws the model faces and mesh edges separately. Each edge is
terminated explicitly, preventing unrelated triangles from being connected.

The public interface and object MAT-file format are unchanged. Existing calls
to `new_object` and `flypath` continue to work, including numeric edge colors
and `'edge','none'`.

The implementation uses MATLAB functionality available in R2014a, the release
against which the original package was tested.

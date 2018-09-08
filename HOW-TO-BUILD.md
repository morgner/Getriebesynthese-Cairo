
Dependencies
============

-- 'c++ 17'
-- 'pthreads'
-- 'gtkmm-3.0'
-- 'cairomm-1.0'


Building the application
========================

```bash
cd <project_root>
mkdir -p build
cd build
cmake ..
cmake --build . --target all -- -j2 (amount of cpu's to use)
```
------------------- OR (e.g.)
```bash
mkdir -p build/CodeBlocks
cd build/CodeBlocks
cmake ../../ -G "CodeBlocks - Unix Makefiles"
codeblocks odb.cbp
```
------------------- OR

use a CMAKE Generator suitable to you environment

Running the application
=======================

./synthese

The application searches for a template in directory ../templates and if found
uses it to generate the file ./synthese.scad which contains the synthesized
Koppelgetriebe.

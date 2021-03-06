#!/usr/bin/python

#  Copyright (C) Vladimir Prus 2003. Permission to copy, use, modify, sell and
#  distribute this software is granted provided this copyright notice appears in
#  all copies. This software is provided "as is" without express or implied
#  warranty, and with no claim as to its suitability for any purpose.

#  Test the 'glob' rule in Jamfile context.
from BoostBuild import Tester, List
import os
import string

# Create a temporary working directory
t = Tester()

t.write("project-root.jam", """ 
""")

t.write("Jamfile", """ 
""")

t.write("d1/a.cpp", """ 
int main() { return 0; }

""")

t.write("d1/Jamfile", """ 
exe a : [ glob *.cpp ] ../d2/d//l ; 
""")

t.write("d2/d/l.cpp", """ 
#if defined(_WIN32)
__declspec(dllexport)
void force_import_lib_creation() {}
#endif
""")

t.write("d2/d/Jamfile", """ 
lib l : [ glob *.cpp ] ; 
""")

t.write("d3/d/Jamfile", """
exe a : [ glob ../*.cpp ] ;
""")
t.write("d3/a.cpp", """
int main()
{
    return 0;
}
""")

t.run_build_system(subdir="d1")
t.expect_addition("d1/bin/$toolset/debug/a.exe")

t.run_build_system(subdir="d3/d")
t.expect_addition("d3/d/bin/$toolset/debug/a.exe")

t.rm("d2/d/bin")
t.run_build_system(subdir="d2/d")
t.expect_addition("d2/d/bin/$toolset/debug/l.dll")

# Test that when 'source-location' is explicitly-specified
# glob works relatively to source location
t.rm("d1")

t.write("d1/src/a.cpp", """ 
int main() { return 0; }

""")

t.write("d1/Jamfile", """
project : source-location src ;
exe a : [ glob *.cpp ] ../d2/d//l ; 
""")

t.run_build_system(subdir="d1")
t.expect_addition("d1/bin/$toolset/debug/a.exe")

# Test that wildcards can include directories
t.rm("d1")

t.write("d1/src/foo/a.cpp", """
void bar();
int main() { bar(); return 0; }

""")

t.write("d1/src/bar/b.cpp", """
void bar() {}

""")


t.write("d1/Jamfile", """
project : source-location src ;
exe a : [ glob foo/*.cpp bar/*.cpp ] ../d2/d//l ; 
""")

t.run_build_system(subdir="d1")
t.expect_addition("d1/bin/$toolset/debug/a.exe")

# Test that 'glob' works with absolute names
t.rm("d1/bin")

# Note that to get current dir, we use bjam's PWD,
# not Python's os.getcwd, because the former will
# always return long path. The latter might return
# short path, and that will confuse path.glob.
t.write("d1/Jamfile", """
project : source-location src ;
local pwd = [ PWD ] ; # Always absolute
exe a : [ glob $(pwd)/src/foo/*.cpp $(pwd)/src/bar/*.cpp ] ../d2/d//l ; 
""")

t.run_build_system(subdir="d1")
t.expect_addition("d1/bin/$toolset/debug/a.exe")


t.cleanup()

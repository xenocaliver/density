# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.29

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/homebrew/Cellar/cmake/3.29.3/bin/cmake

# The command to remove a file.
RM = /opt/homebrew/Cellar/cmake/3.29.3/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/jason/workspace/density

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/jason/workspace/density/build

# Include any dependencies generated for this target.
include CMakeFiles/density.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/density.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/density.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/density.dir/flags.make

density_autogen/timestamp: /opt/homebrew/share/qt/libexec/moc
density_autogen/timestamp: /opt/homebrew/share/qt/libexec/uic
density_autogen/timestamp: CMakeFiles/density.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --blue --bold --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Automatic MOC and UIC for target density"
	/opt/homebrew/Cellar/cmake/3.29.3/bin/cmake -E cmake_autogen /Users/jason/workspace/density/build/CMakeFiles/density_autogen.dir/AutogenInfo.json Debug
	/opt/homebrew/Cellar/cmake/3.29.3/bin/cmake -E touch /Users/jason/workspace/density/build/density_autogen/timestamp

CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o: CMakeFiles/density.dir/flags.make
CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o: density_autogen/mocs_compilation.cpp
CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o: CMakeFiles/density.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o -MF CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o.d -o CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o -c /Users/jason/workspace/density/build/density_autogen/mocs_compilation.cpp

CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jason/workspace/density/build/density_autogen/mocs_compilation.cpp > CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.i

CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jason/workspace/density/build/density_autogen/mocs_compilation.cpp -o CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.s

CMakeFiles/density.dir/density.cpp.o: CMakeFiles/density.dir/flags.make
CMakeFiles/density.dir/density.cpp.o: /Users/jason/workspace/density/density.cpp
CMakeFiles/density.dir/density.cpp.o: CMakeFiles/density.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/density.dir/density.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/density.dir/density.cpp.o -MF CMakeFiles/density.dir/density.cpp.o.d -o CMakeFiles/density.dir/density.cpp.o -c /Users/jason/workspace/density/density.cpp

CMakeFiles/density.dir/density.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/density.dir/density.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jason/workspace/density/density.cpp > CMakeFiles/density.dir/density.cpp.i

CMakeFiles/density.dir/density.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/density.dir/density.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jason/workspace/density/density.cpp -o CMakeFiles/density.dir/density.cpp.s

CMakeFiles/density.dir/plot.cpp.o: CMakeFiles/density.dir/flags.make
CMakeFiles/density.dir/plot.cpp.o: /Users/jason/workspace/density/plot.cpp
CMakeFiles/density.dir/plot.cpp.o: CMakeFiles/density.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/density.dir/plot.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/density.dir/plot.cpp.o -MF CMakeFiles/density.dir/plot.cpp.o.d -o CMakeFiles/density.dir/plot.cpp.o -c /Users/jason/workspace/density/plot.cpp

CMakeFiles/density.dir/plot.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/density.dir/plot.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jason/workspace/density/plot.cpp > CMakeFiles/density.dir/plot.cpp.i

CMakeFiles/density.dir/plot.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/density.dir/plot.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jason/workspace/density/plot.cpp -o CMakeFiles/density.dir/plot.cpp.s

CMakeFiles/density.dir/load_json.cpp.o: CMakeFiles/density.dir/flags.make
CMakeFiles/density.dir/load_json.cpp.o: /Users/jason/workspace/density/load_json.cpp
CMakeFiles/density.dir/load_json.cpp.o: CMakeFiles/density.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/density.dir/load_json.cpp.o"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/density.dir/load_json.cpp.o -MF CMakeFiles/density.dir/load_json.cpp.o.d -o CMakeFiles/density.dir/load_json.cpp.o -c /Users/jason/workspace/density/load_json.cpp

CMakeFiles/density.dir/load_json.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Preprocessing CXX source to CMakeFiles/density.dir/load_json.cpp.i"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /Users/jason/workspace/density/load_json.cpp > CMakeFiles/density.dir/load_json.cpp.i

CMakeFiles/density.dir/load_json.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green "Compiling CXX source to assembly CMakeFiles/density.dir/load_json.cpp.s"
	/Library/Developer/CommandLineTools/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /Users/jason/workspace/density/load_json.cpp -o CMakeFiles/density.dir/load_json.cpp.s

# Object files for target density
density_OBJECTS = \
"CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o" \
"CMakeFiles/density.dir/density.cpp.o" \
"CMakeFiles/density.dir/plot.cpp.o" \
"CMakeFiles/density.dir/load_json.cpp.o"

# External object files for target density
density_EXTERNAL_OBJECTS =

density: CMakeFiles/density.dir/density_autogen/mocs_compilation.cpp.o
density: CMakeFiles/density.dir/density.cpp.o
density: CMakeFiles/density.dir/plot.cpp.o
density: CMakeFiles/density.dir/load_json.cpp.o
density: CMakeFiles/density.dir/build.make
density: /opt/homebrew/lib/QtWidgets.framework/Versions/A/QtWidgets
density: /opt/homebrew/lib/qwt.framework/qwt
density: lib/libfftw3.a
density: lib/libfftw3_threads.a
density: /opt/homebrew/lib/QtGui.framework/Versions/A/QtGui
density: /opt/homebrew/lib/QtCore.framework/Versions/A/QtCore
density: CMakeFiles/density.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color "--switch=$(COLOR)" --green --bold --progress-dir=/Users/jason/workspace/density/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable density"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/density.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/density.dir/build: density
.PHONY : CMakeFiles/density.dir/build

CMakeFiles/density.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/density.dir/cmake_clean.cmake
.PHONY : CMakeFiles/density.dir/clean

CMakeFiles/density.dir/depend: density_autogen/timestamp
	cd /Users/jason/workspace/density/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/jason/workspace/density /Users/jason/workspace/density /Users/jason/workspace/density/build /Users/jason/workspace/density/build /Users/jason/workspace/density/build/CMakeFiles/density.dir/DependInfo.cmake "--color=$(COLOR)"
.PHONY : CMakeFiles/density.dir/depend

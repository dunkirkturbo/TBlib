# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.20

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "C:\Program Files\JetBrains\CLion 2021.2.2\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "C:\Program Files\JetBrains\CLion 2021.2.2\bin\cmake\win\bin\cmake.exe" -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = C:\Users\lenovo\Desktop\TBlib

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = C:\Users\lenovo\Desktop\TBlib\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/TBlib.dir/depend.make
# Include the progress variables for this target.
include CMakeFiles/TBlib.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/TBlib.dir/flags.make

CMakeFiles/TBlib.dir/main.cpp.obj: CMakeFiles/TBlib.dir/flags.make
CMakeFiles/TBlib.dir/main.cpp.obj: CMakeFiles/TBlib.dir/includes_CXX.rsp
CMakeFiles/TBlib.dir/main.cpp.obj: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\lenovo\Desktop\TBlib\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/TBlib.dir/main.cpp.obj"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\TBlib.dir\main.cpp.obj -c C:\Users\lenovo\Desktop\TBlib\main.cpp

CMakeFiles/TBlib.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TBlib.dir/main.cpp.i"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\lenovo\Desktop\TBlib\main.cpp > CMakeFiles\TBlib.dir\main.cpp.i

CMakeFiles/TBlib.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TBlib.dir/main.cpp.s"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\lenovo\Desktop\TBlib\main.cpp -o CMakeFiles\TBlib.dir\main.cpp.s

CMakeFiles/TBlib.dir/utils/bench.cpp.obj: CMakeFiles/TBlib.dir/flags.make
CMakeFiles/TBlib.dir/utils/bench.cpp.obj: CMakeFiles/TBlib.dir/includes_CXX.rsp
CMakeFiles/TBlib.dir/utils/bench.cpp.obj: ../utils/bench.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=C:\Users\lenovo\Desktop\TBlib\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/TBlib.dir/utils/bench.cpp.obj"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\TBlib.dir\utils\bench.cpp.obj -c C:\Users\lenovo\Desktop\TBlib\utils\bench.cpp

CMakeFiles/TBlib.dir/utils/bench.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/TBlib.dir/utils/bench.cpp.i"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E C:\Users\lenovo\Desktop\TBlib\utils\bench.cpp > CMakeFiles\TBlib.dir\utils\bench.cpp.i

CMakeFiles/TBlib.dir/utils/bench.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/TBlib.dir/utils/bench.cpp.s"
	C:\PROGRA~1\mingw64\bin\G__~1.EXE $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S C:\Users\lenovo\Desktop\TBlib\utils\bench.cpp -o CMakeFiles\TBlib.dir\utils\bench.cpp.s

# Object files for target TBlib
TBlib_OBJECTS = \
"CMakeFiles/TBlib.dir/main.cpp.obj" \
"CMakeFiles/TBlib.dir/utils/bench.cpp.obj"

# External object files for target TBlib
TBlib_EXTERNAL_OBJECTS =

TBlib.exe: CMakeFiles/TBlib.dir/main.cpp.obj
TBlib.exe: CMakeFiles/TBlib.dir/utils/bench.cpp.obj
TBlib.exe: CMakeFiles/TBlib.dir/build.make
TBlib.exe: libNTL.a
TBlib.exe: libTBlib_extra.dll.a
TBlib.exe: CMakeFiles/TBlib.dir/linklibs.rsp
TBlib.exe: CMakeFiles/TBlib.dir/objects1.rsp
TBlib.exe: CMakeFiles/TBlib.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=C:\Users\lenovo\Desktop\TBlib\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CXX executable TBlib.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\TBlib.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/TBlib.dir/build: TBlib.exe
.PHONY : CMakeFiles/TBlib.dir/build

CMakeFiles/TBlib.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\TBlib.dir\cmake_clean.cmake
.PHONY : CMakeFiles/TBlib.dir/clean

CMakeFiles/TBlib.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" C:\Users\lenovo\Desktop\TBlib C:\Users\lenovo\Desktop\TBlib C:\Users\lenovo\Desktop\TBlib\cmake-build-debug C:\Users\lenovo\Desktop\TBlib\cmake-build-debug C:\Users\lenovo\Desktop\TBlib\cmake-build-debug\CMakeFiles\TBlib.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/TBlib.dir/depend


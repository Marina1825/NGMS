# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

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
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/marina/4_curs/NGMS

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/marina/4_curs/NGMS/build

# Include any dependencies generated for this target.
include CMakeFiles/tcp_proxy.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/tcp_proxy.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/tcp_proxy.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/tcp_proxy.dir/flags.make

CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o: CMakeFiles/tcp_proxy.dir/flags.make
CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o: ../libmzmq/tcp_proxy.cpp
CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o: CMakeFiles/tcp_proxy.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/marina/4_curs/NGMS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o -MF CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o.d -o CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o -c /home/marina/4_curs/NGMS/libmzmq/tcp_proxy.cpp

CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/marina/4_curs/NGMS/libmzmq/tcp_proxy.cpp > CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.i

CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/marina/4_curs/NGMS/libmzmq/tcp_proxy.cpp -o CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.s

# Object files for target tcp_proxy
tcp_proxy_OBJECTS = \
"CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o"

# External object files for target tcp_proxy
tcp_proxy_EXTERNAL_OBJECTS =

tcp_proxy: CMakeFiles/tcp_proxy.dir/libmzmq/tcp_proxy.cpp.o
tcp_proxy: CMakeFiles/tcp_proxy.dir/build.make
tcp_proxy: CMakeFiles/tcp_proxy.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/marina/4_curs/NGMS/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable tcp_proxy"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/tcp_proxy.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/tcp_proxy.dir/build: tcp_proxy
.PHONY : CMakeFiles/tcp_proxy.dir/build

CMakeFiles/tcp_proxy.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/tcp_proxy.dir/cmake_clean.cmake
.PHONY : CMakeFiles/tcp_proxy.dir/clean

CMakeFiles/tcp_proxy.dir/depend:
	cd /home/marina/4_curs/NGMS/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/marina/4_curs/NGMS /home/marina/4_curs/NGMS /home/marina/4_curs/NGMS/build /home/marina/4_curs/NGMS/build /home/marina/4_curs/NGMS/build/CMakeFiles/tcp_proxy.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/tcp_proxy.dir/depend


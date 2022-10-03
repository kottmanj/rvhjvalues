cmake_minimum_required (VERSION 3.0)
project(h2j)

enable_language(CXX)
include(GNUInstallDirs)

# need MADNESS
find_package(MADNESS 0.10.1 CONFIG REQUIRED QUIET COMPONENTS chem HINTS ${MAD_ROOT_DIR})

# the main executable
add_executable (h2j h2j.cpp)
#target_include_directories(h2j PUBLIC "/home/jsk/anaconda3/envs/tq-1.8.1/")
target_link_libraries(h2j LINK_PUBLIC MADchem)



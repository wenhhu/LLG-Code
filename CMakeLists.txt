cmake_minimum_required (VERSION 3.5)
project (Magnonic_superlattice)

include_directories (eigen3)
include_directories (include)
include_directories (“/usr/local/include/gsl”)

set ( PROJECT_LINK_LIBS libgsl.a )

file(GLOB SOURCES "src/*.cpp")

add_executable(Magnon ${SOURCES})
target_link_libraries(Magnon ${PROJECT_LINK_LIBS})

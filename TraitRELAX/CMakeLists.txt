# CMake script for traitrelax main programs
# Author:
#   Keren Halabi
#   Francois Gindraud (2017)
# Created: 02/08/2020

# Add all programs.
# Any .cpp file in TraitRELAX/ is considered to be a program.
# It will be compiled as a standalone program (must contain a main()).

file (GLOB program_cpp_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)
foreach (program_cpp_file ${program_cpp_files})
  # Add each test (named as the filename without extension)
  get_filename_component (program_name ${program_cpp_file} NAME_WE)
  add_executable (${program_name} ${program_cpp_file})
  target_link_libraries (${program_name} ${PROJECT_NAME}-shared)
  set_target_properties (${program_name} PROPERTIES POSITION_INDEPENDENT_CODE TRUE)
  add_test (
    NAME ${program_name}
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMMAND ${program_name}
    )
endforeach (program_cpp_file)

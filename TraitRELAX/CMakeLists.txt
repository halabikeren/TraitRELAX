# CMake script for TraitRELAX program, based on the script from TestNH (https://github.com/BioPP/testnh)

# Executables of TraitRELAX
add_executable (traitrelax TraitRELAX.cpp)

set (traitrelax-targets
  traitrelax
  )

foreach (target ${traitrelax-targets})
  # Link (static or shared)
  if (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_STATIC})
    set_target_properties (${target} PROPERTIES LINK_SEARCH_END_STATIC TRUE)
  else (BUILD_STATIC)
    target_link_libraries (${target} ${BPP_LIBS_SHARED})
  endif (BUILD_STATIC)
endforeach (target)

install (TARGETS ${testnh-targets} DESTINATION ${CMAKE_INSTALL_BINDIR})

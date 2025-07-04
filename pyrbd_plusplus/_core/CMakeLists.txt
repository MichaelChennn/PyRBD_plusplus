# set source files
set(SRC
    mcs.cpp
    pathset.cpp
    sdp.cpp
    utils.cpp
    bindings.cpp
)

# find pybind11
include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG        v2.13.6
)
FetchContent_MakeAvailable(pybind11)

# find OpenMP
find_package(OpenMP REQUIRED)

# create the library
add_library(pyrbd_plusplus MODULE ${SRC})


# link the pybind11 library
target_link_libraries(pyrbd_plusplus PRIVATE pybind11::module)

# set complie options
target_compile_options(pyrbd_plusplus PRIVATE
    -Wall -Wextra -O3 -march=native
)

# Link OpenMP if found
if(OpenMP_CXX_FOUND)
    target_link_libraries(pyrbd_plusplus PUBLIC OpenMP::OpenMP_CXX)
endif()

# set the include directories
target_include_directories(pyrbd_plusplus PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

# set the executable properties
set_target_properties(pyrbd_plusplus PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    OUTPUT_NAME "pyrbd_plusplus"
    PREFIX ""
    SUFFIX ".so"
)

# clean target
add_custom_target(distclean
    COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
    COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt
    COMMAND ${CMAKE_COMMAND} -E remove cmake_install.cmake
    COMMAND ${CMAKE_COMMAND} -E remove Makefile
    COMMAND ${CMAKE_COMMAND} -E remove -f ${CMAKE_CURRENT_SOURCE_DIR}/*.so
    COMMAND ${CMAKE_COMMAND} -E remove_directory _deps
    COMMENT "Removing all build files"
)
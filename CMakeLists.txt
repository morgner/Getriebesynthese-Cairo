cmake_minimum_required(VERSION 3.9)
project(synthese
  LANGUAGES CXX
  VERSION 0.1.20
  )
 
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra -Werror -pedantic-errors")

set(CMAKE_INTERPROCEDURAL_OPTIMIZATION ON)
set(CMAKE_CXX_EXTENSIONS OFF)
#set(CMAKE_CXX_STANDARD "2a")
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

find_package("Threads")
find_package("PkgConfig")

pkg_check_modules(GTK
    REQUIRED
    "gtkmm-3.0") #gtkmm-3.0
#    "gtk+-3.0") #gtkmm-3.0
pkg_check_modules(CAIRO
    REQUIRED
    "cairomm-1.0") #cairomm-1.0
#    "cairo") #cairomm-1.0

include_directories(SYSTEM
    ${CAIRO_INCLUDE_DIRS}
    ${GTK_INCLUDE_DIRS}
)
include_directories("include")

add_subdirectory("source")
add_subdirectory("include")

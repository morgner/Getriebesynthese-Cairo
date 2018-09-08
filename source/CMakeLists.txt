#cmake_minimum_required(VERSION 3.6)
#project(synthese.cairo)
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/cmake")
#set(CMAKE_CXX_STANDARD 17)
#
#set(SOURCE_FILES main.cpp)
#add_executable(synthese.cairo ${SOURCE_FILES})
#
##find_package(ImageMagick COMPONENTS Magick++)
##include_directories(${ImageMagick_INCLUDE_DIRS})
##target_link_libraries(HelloOpenGL ${ImageMagick_LIBRARIES})
#
#find_package(Cairo)
#include_directories(${Cairo_INCLUDE_DIRS})
#target_link_libraries(synthese.cairo ${Cairo_LIBRARIES})

add_executable(${PROJECT_NAME}
    "main.cpp"
    "canvas.cpp"
)
target_link_libraries(${PROJECT_NAME}
    ${CMAKE_THREAD_LIBS_INIT}
    ${CAIRO_LIBRARIES}
    ${GTK_LIBRARIES}
)

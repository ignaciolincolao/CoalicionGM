# CMakeList.txt: proyecto de CMake para CGM, incluya el origen y defina
# la lógica específica del proyecto aquí.
#
cmake_minimum_required (VERSION 3.8)

project ("CGM")

include_directories(
	${PROJECT_SOURCE_DIR}/include
)

file(GLOB all_SRCS
	"${PROJECT_SOURCE_DIR}/src/*.cpp"
	"${PROJECT_SOURCE_DIR}/src/*.c"
	"${PROJECT_SOURCE_DIR}/src/*.h"
)
file(GLOB all_DATA
	"${PROJECT_SOURCE_DIR}/data/*.json"
)
# Agregue un origen al ejecutable de este proyecto.
add_executable (CGM ${all_SRCS} ${all_DATA})
set_target_properties(CGM PROPERTIES LINKER_LANGUAGE CXX)

file(COPY ${all_DATA} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#configure_file(${CMAKE_CURRENT_SOURCE_DIR}/data/votacion.json
	#${CMAKE_CURRENT_BINARY_DIR} COPYONLY)
# TODO: Agregue pruebas y destinos de instalación si es necesario.

include(FetchContent)

FetchContent_Declare(json URL https://github.com/nlohmann/json/releases/download/v3.10.5/json.tar.xz)
FetchContent_MakeAvailable(json)

target_link_libraries(CGM PRIVATE nlohmann_json::nlohmann_json)


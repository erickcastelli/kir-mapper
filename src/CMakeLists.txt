
#BOOST_ROOT=path_to_conda_env ZLIB_ROOT=path_to_conda_env cmake ..

cmake_minimum_required (VERSION 3.5)
#cmake_policy (SET CMP0144 NEW)
cmake_policy (SET CMP0074 NEW)

SET (CMAKE_INSTALL_PREFIX "/usr/local/bin")
SET(BUILD_SHARED_LIBS OFF)
SET(CMAKE_EXE_LINKER_FLAGS "-O3 -std=c++17 -lboost_iostreams -lz -L /usr/local/lib -I /usr/local/include")
#set(Boost_DEBUG 1)

project (kir-mapper)

SET(Boost_USE_STATIC_LIBS   ON)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME ON)
set(Boost_NO_BOOST_CMAKE     ON)

message("")
find_package(ZLIB REQUIRED)
message(ZLIB_INCLUDE   - [${ZLIB_INCLUDE_DIRS}])
message(ZLIB_LIBRARIES - [${ZLIB_LIBRARIES}])

message("")
find_package( Boost 1.70 COMPONENTS program_options regex system iostreams filesystem REQUIRED)
message(Boost_LIBRARIES   - [${Boost_LIBRARY_DIRS}])
message(Boost_INCLUDE   - [${Boost_INCLUDE_DIRS}])

message("")
set (CMAKE_CXX_STANDARD 17)
add_executable(kir-mapper external.hpp setup.hpp setup.cpp combine.cpp combine.hpp exome_bias.cpp exome_bias.hpp functions.cpp functions.hpp main.cpp genotype.cpp genotype.hpp map_dna.cpp map_dna.hpp ThreadPool.hpp haplotypes.cpp haplotypes.hpp ncopy.cpp ncopy.hpp preselect.cpp preselect.hpp)
target_link_libraries(kir-mapper pthread Boost::program_options Boost::regex Boost::system Boost::iostreams Boost::filesystem ${ZLIB_LIBRARIES})
#target_link_libraries(kir-mapper pthread ${BOOST_LIBRARIES} ${ZLIB_LIBRARIES})

cmake_minimum_required(VERSION 3.5)
project(TranSurVeyor)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -pthread")

include_directories(htslib-1.11/include)

set(READS_CATEGORIZER_SRC reads_categorizer.cpp sam_utils.h config.h libs/cptl_stl.h)
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})

set(CLIP_CONSENSUS_BUILDER_SRC libs/ssw.c libs/ssw_cpp.cpp clip_consensus_builder.cpp sam_utils.h config.h libs/cptl_stl.h sw_utils.h cluster.h)
add_executable(clip_consensus_builder ${CLIP_CONSENSUS_BUILDER_SRC})

set(SC_CATEGORIZER_SRC sc_categorizer.cpp sam_utils.h config.h libs/cptl_stl.h cluster.h)
add_executable(sc_categorizer ${SC_CATEGORIZER_SRC})

set(DC_REMAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp dc_remapper.cpp sam_utils.h config.h libs/cptl_stl.h cluster.h)
add_executable(dc_remapper ${DC_REMAPPER_SRC})

set(ADD_FILTERING_INFO_SRC add_filtering_info.cpp sam_utils.h config.h libs/cptl_stl.h cluster.h ks-test.h)
add_executable(add_filtering_info ${ADD_FILTERING_INFO_SRC})

set(FILTER_SRC filter.cpp cluster.h)
add_executable(filter ${FILTER_SRC})

find_library(HTS_LIB hts PATHS htslib-1.11/lib NO_DEFAULT_PATH)
target_link_libraries(reads_categorizer "${HTS_LIB}")
target_link_libraries(clip_consensus_builder "${HTS_LIB}")
target_link_libraries(sc_categorizer "${HTS_LIB}")
target_link_libraries(dc_remapper "${HTS_LIB}")
target_link_libraries(add_filtering_info "${HTS_LIB}")
target_link_libraries(filter "${HTS_LIB}")

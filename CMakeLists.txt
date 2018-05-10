cmake_minimum_required(VERSION 3.5)
project(TranSurVeyor)

set(CMAKE_CXX_STANDARD 11)

include_directories(htslib-1.7/include)
link_directories(htslib-1.7/lib)

set(HEADERS_SRC sam_utils.h config.h libs/cptl_stl.h sw_utils.h cluster.h)

set(READS_CATEGORIZER_SRC reads_categorizer.cpp ${HEADERS_SRC})
add_executable(reads_categorizer ${READS_CATEGORIZER_SRC})
target_link_libraries(reads_categorizer hts)

set(CLIP_CONSENSUS_BUILDER_SRC libs/ssw.c libs/ssw_cpp.cpp clip_consensus_builder.cpp ${HEADERS_SRC})
add_executable(clip_consensus_builder ${CLIP_CONSENSUS_BUILDER_SRC})
target_link_libraries(clip_consensus_builder hts)

set(SC_CATEGORIZER_SRC sc_categorizer.cpp ${HEADERS_SRC})
add_executable(sc_categorizer ${SC_CATEGORIZER_SRC})
target_link_libraries(sc_categorizer hts)

set(DC_REMAPPER_SRC libs/ssw.c libs/ssw_cpp.cpp dc_remapper.cpp ${HEADERS_SRC})
add_executable(dc_remapper ${DC_REMAPPER_SRC})
target_link_libraries(dc_remapper hts)

set(CLUSTERER_SRC clusterer.cpp ${HEADERS_SRC})
add_executable(clusterer ${CLUSTERER_SRC})
target_link_libraries(clusterer hts)

set(ADD_FILTERING_INFO_SRC add_filtering_info.cpp ${HEADERS_SRC})
add_executable(add_filtering_info ${ADD_FILTERING_INFO_SRC})
target_link_libraries(add_filtering_info hts)

set(FILTER_SRC filter.cpp ${HEADERS_SRC})
add_executable(filter ${FILTER_SRC})
target_link_libraries(filter hts)

#find_library(HTS_LIB hts)
#target_link_libraries(reads_categorizer "${HTS_LIB}")
#target_link_libraries(clip_consensus_builder "${HTS_LIB}")
#target_link_libraries(sc_categorizer "${HTS_LIB}")
#target_link_libraries(dc_remapper "${HTS_LIB}")
#target_link_libraries(clusterer "${HTS_LIB}")
##target_link_libraries(add_filtering_info "${HTS_LIB}")
#target_link_libraries(filter "${HTS_LIB}")

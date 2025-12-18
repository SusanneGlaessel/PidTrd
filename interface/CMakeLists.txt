#include(AnalysisTree)

set(SOURCES
  PidTrdMcInput.cpp
  PidTrdRunGetter.cpp
  PidTrdFiller.cpp
  )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
add_library(PidTrdMan SHARED ${SOURCES} G__PidTrdMan.cxx)

#message(${PROJECT_INCLUDE_DIRECTORIES})

#include_directories(${PROJECT_INCLUDE_DIRECTORIES} ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src)
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${AnalysisTree_INCLUDE_DIR})

ROOT_GENERATE_DICTIONARY(G__PidTrdMan ${HEADERS} LINKDEF PidTrdManLinkDef.h)
target_link_libraries(PidTrdMan ${ROOT_LIBRARIES} AnalysisTreeBase AnalysisTreeInfra PidTrd)

add_executable(main main.cpp)
target_link_libraries(main PidTrd PidTrdMan AnalysisTreeBase AnalysisTreeInfra ${ROOT_LIBRARIES})
target_include_directories(main PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/interface ${PROJECT_INCLUDE_DIRECTORIES})

install(TARGETS main RUNTIME DESTINATION bin)

add_executable(create_mcinput create_mcinput.cpp)
target_link_libraries(create_mcinput PidTrd PidTrdMan AnalysisTreeBase AnalysisTreeInfra ${ROOT_LIBRARIES})
target_include_directories(create_mcinput PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/interface ${PROJECT_INCLUDE_DIRECTORIES})

install(TARGETS create_mcinput RUNTIME DESTINATION bin)

add_executable(create_getter create_getter.cpp)
target_link_libraries(create_getter PidTrd PidTrdMan ${ROOT_LIBRARIES})
target_include_directories(create_getter PUBLIC ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/interface ${PROJECT_INCLUDE_DIRECTORIES})

install(TARGETS create_getter RUNTIME DESTINATION bin)

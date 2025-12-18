set(SOURCES
        ParticleProb.cpp	
	TrdContainer.cpp
        Getter.cpp
        )

string(REPLACE ".cpp" ".hpp" HEADERS "${SOURCES}")
list(APPEND HEADERS Constants.hpp)

set(DICT_FILE_NAME G__${PROJECT_NAME})
set(PCM_FILE_NAME lib${PROJECT_NAME})

include_directories(${CMAKE_CURRENT_SOURCE_DIR})
ROOT_GENERATE_DICTIONARY(${DICT_FILE_NAME} ${HEADERS} LINKDEF PidTrdLinkDef.h)
add_library(PidTrd SHARED ${SOURCES} ${HEADERS} ${DICT_FILE_NAME})
foreach (_HEADER ${HEADERS})
    get_filename_component(_HEADER_FILENAME ${_HEADER} NAME)
    get_filename_component(_HEADER_PATH ${_HEADER} REALPATH)
    set(_INCLUDE_BINARY_DIR ${CMAKE_BINARY_DIR}/include/pidtrd)
    set(_HEADER_SYMLINK_PATH ${_INCLUDE_BINARY_DIR}/${_HEADER_FILENAME})
    add_custom_command(TARGET PidTrd PRE_BUILD
            COMMAND ${CMAKE_COMMAND} ARGS -E make_directory ${_INCLUDE_BINARY_DIR}
            COMMAND ${CMAKE_COMMAND} ARGS -E create_symlink ${_HEADER_PATH} ${_HEADER_SYMLINK_PATH}
            )
endforeach ()

target_link_libraries(PidTrd ${ROOT_LIBRARIES})
target_include_directories(PidTrd
        PRIVATE ${CMAKE_SOURCE_DIR}/src
        PUBLIC
        $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include/pidtrd>
        $<INSTALL_INTERFACE:include/pidtrd>
        )

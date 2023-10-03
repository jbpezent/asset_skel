include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/extensions)



set(ASSET_INCLUDES
  ${CMAKE_CURRENT_SOURCE_DIR}
  PARENT_SCOPE
)

#-------------------------------------------------------------------------------
## PCH
add_library(pch STATIC pch.h pch.cpp  )
target_precompile_headers(pch PRIVATE pch.h)
target_include_directories(pch PRIVATE ${INCLUDE_DIRS})
target_compile_options(pch PUBLIC ${COMPILE_FLAGS})



#---------------------------------------------------------------------------------
## Static Library
add_library(staticlib 
STATIC 
staticlib.h staticlib.cpp
TestMatrices.h TestMatrices.cpp
)

target_include_directories(staticlib PRIVATE ${INCLUDE_DIRS})
target_compile_options(staticlib PUBLIC ${COMPILE_FLAGS})
target_precompile_headers(staticlib REUSE_FROM pch)

add_executable(testexec main.cpp)
target_compile_options(testexec PUBLIC  ${COMPILE_FLAGS})
target_include_directories(testexec PRIVATE ${INCLUDE_DIRS})
target_link_libraries(testexec PUBLIC staticlib)
target_link_libraries(testexec PRIVATE ${LINK_LIBS})
target_precompile_headers(testexec REUSE_FROM pch)
target_link_libraries(staticlib PRIVATE fmt::fmt-header-only)



#-------------------------------------------------------------------------------
## ASSET Python Module
pybind11_add_module(_asset_skel main.cpp)
target_compile_options(_asset_skel PUBLIC  ${COMPILE_FLAGS})
target_link_libraries(_asset_skel PUBLIC staticlib)
target_include_directories(_asset_skel PRIVATE ${INCLUDE_DIRS})
target_link_libraries(_asset_skel PRIVATE ${LINK_LIBS})

if(UNIX)
    set_target_properties(_asset_skel PROPERTIES LINK_FLAGS ${LINKER_FLAGS})
endif()

add_custom_target(pyassetsrc ALL DEPENDS ${CMAKE_SOURCE_DIR}/asset_skel)


if(PYTHON_LOCAL_INSTALL_DIR)
    add_custom_command(TARGET _asset_skel POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_asset_skel> ${PYTHON_LOCAL_INSTALL_DIR}/
    COMMENT "Copied _asset_skel.pyd to ${PYTHON_LOCAL_INSTALL_DIR}"
    )
    message(STATUS "Cleaning pyasset install directory")

    add_custom_command(TARGET pyassetsrc POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E rm -rf ${PYTHON_LOCAL_INSTALL_DIR}/asset_skel
                    COMMENT "Cleaning asset_skel install directory")

    add_custom_command(TARGET pyassetsrc POST_BUILD
                    COMMAND ${CMAKE_COMMAND} -E copy_directory  ${CMAKE_SOURCE_DIR}/asset_skel
                                                ${PYTHON_LOCAL_INSTALL_DIR}/asset_skel
                    DEPENDS ${CMAKE_SOURCE_DIR}/asset_skel
                    COMMENT "Copied asset_skel to ${PYTHON_LOCAL_INSTALL_DIR}")

endif()


    




    



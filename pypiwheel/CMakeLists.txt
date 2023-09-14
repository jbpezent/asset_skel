
set(CMAKE_VERBOSE_MAKEFILE TRUE)

add_custom_target(pypiwheel  ALL DEPENDS ${CMAKE_SOURCE_DIR}/asset_skel)
add_dependencies(pypiwheel _asset_skel)

set(WHEEL_BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/asset_skel)

set(ASSET_BINARY_NAME "_asset_skel${PYTHON_MODULE_EXTENSION}")


configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/setup.py.in"
    "${CMAKE_CURRENT_BINARY_DIR}/setup.py"
    @ONLY)


add_custom_command(TARGET pypiwheel 
                    COMMAND ${CMAKE_COMMAND} -E rm -rf ${WHEEL_BUILD_DIR}/
                    COMMENT "Cleaning asset_asrl install directory")


add_custom_command(TARGET pypiwheel 
                    COMMAND ${CMAKE_COMMAND} -E copy_directory  ${CMAKE_SOURCE_DIR}/asset_skel
                                                ${WHEEL_BUILD_DIR}/asset_skel
                    DEPENDS ${CMAKE_SOURCE_DIR}/asset_skel
                    COMMENT "Copied asset_asrl to ${WHEEL_BUILD_DIR}")

add_custom_command(TARGET pypiwheel 
                    COMMAND ${CMAKE_COMMAND} -E copy  ${CMAKE_CURRENT_BINARY_DIR}/setup.py ${WHEEL_BUILD_DIR}
                    
                    COMMENT "Copied setup.py to ${WHEEL_BUILD_DIR}")


add_custom_command(TARGET pypiwheel 
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:_asset_skel> ${WHEEL_BUILD_DIR}/
    COMMENT "Copied asset.pyd to ${WHEEL_BUILD_DIR}"
    )




add_custom_command(TARGET pypiwheel 
                    COMMAND
                   ${PYTHON_EXECUTABLE} ${WHEEL_BUILD_DIR}/setup.py bdist_wheel 
                   WORKING_DIRECTORY ${WHEEL_BUILD_DIR}
                   COMMENT "Generating Binary Wheel")


message(STATUS ${ASSET_BINARY_NAME})



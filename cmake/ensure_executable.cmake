if(DEFINED TARGET_FILES AND NOT TARGET_FILES STREQUAL "")
    foreach(_target_file IN LISTS TARGET_FILES)
        if(EXISTS "${_target_file}")
            file(CHMOD "${_target_file}"
                 PERMISSIONS
                 OWNER_READ OWNER_WRITE OWNER_EXECUTE
                 GROUP_READ GROUP_EXECUTE
                 WORLD_READ WORLD_EXECUTE)
        endif()
    endforeach()
elseif(DEFINED TARGET_FILE AND NOT TARGET_FILE STREQUAL "")
    file(CHMOD "${TARGET_FILE}"
         PERMISSIONS
         OWNER_READ OWNER_WRITE OWNER_EXECUTE
         GROUP_READ GROUP_EXECUTE
         WORLD_READ WORLD_EXECUTE)
else()
    message(FATAL_ERROR "TARGET_FILE or TARGET_FILES is required")
endif()

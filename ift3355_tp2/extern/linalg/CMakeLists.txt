add_library(linalg "")
set_target_properties(linalg PROPERTIES LINKER_LANGUAGE CXX)


target_sources(linalg
                PUBLIC
        linalg.h
                )

target_include_directories(linalg
                            PUBLIC
        .
                            )
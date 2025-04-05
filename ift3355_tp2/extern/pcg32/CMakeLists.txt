add_library(pcg32 "")
set_target_properties(pcg32 PROPERTIES LINKER_LANGUAGE CXX)


target_sources(pcg32
                PUBLIC
        pcg32.h
                )

target_include_directories(pcg32
                            PUBLIC
        .
                        )
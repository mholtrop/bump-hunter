
add_library(HpsFitResult HpsFitResult.cxx
  ${ROOT_LIBRARIES}
)

add_library(BumpHunter BumpHunter.cxx)
target_link_libraries(BumpHunter
    ${ROOT_LIBRARIES}
    HpsFitResult
    FitPrinter
)

add_executable(bump_hunter bump_hunter.cxx)
target_link_libraries(bump_hunter
    ${ROOT_LIBRARIES}
    BumpHunter
    FlatTupleMaker
    HpsFitResult
)

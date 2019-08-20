
add_library(RootFileReader RootFileReader.cxx)
target_link_libraries(RootFileReader 
   ${ROOT_LIBRARIES}
)

add_library(FlatTupleMaker FlatTupleMaker.cxx)
target_link_libraries(FlatTupleMaker
   ${ROOT_LIBRARIES}
)

add_library(FitPrinter FitPrinter.cxx)
target_link_libraries(FitPrinter
   ${ROOT_LIBRARIES}
)

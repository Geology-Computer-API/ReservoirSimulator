########### next target ###############

file(GLOB HEADERS *.h)

file(GLOB SOURCE *.cpp)

add_executable(MalhasHdiv MalhasHdiv.cpp ${lib_files} ${SOURCE} ${HEADERS})



target_include_directories(MalhasHdiv PUBLIC ${lib_dir})

target_link_libraries(MalhasHdiv pz ${OpenCV_LIBS})
file (COPY PlotGenerator.nb/ DESTINATION Debug/PlotGenerator.nb)
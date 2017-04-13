cmake_minimum_required(VERSION 2.6)

PROJECT(CorrectAttenuation_SV)

FIND_PACKAGE(ITK)
IF(ITK_FOUND)
INCLUDE(${ITK_USE_FILE})
ELSE(ITK_FOUND)
MESSAGE(FATAL_ERROR
"ITK not found. Please set ITK_DIR.")
ENDIF(ITK_FOUND)

#find_package(ItkVtkGlue REQUIRED)
#include(${ItkVtkGlue_USE_FILE})


FIND_PACKAGE ( VTK)
IF (   VTK_FOUND   )
         INCLUDE(  ${VTK_USE_FILE}  )
ELSE(VTK_FOUND)
MESSAGE(FATAL_ERROR
"VTK not found. Please set VTK_DIR.")
ENDIF(   VTK_FOUND   )


# Build a static executable
if (STATIC_BUILD)
 SET(CMAKE_FIND_LIBRARY_SUFFIXES ".a")
 SET(BUILD_SHARED_LIBRARIES OFF)
 SET(CMAKE_EXE_LINKER_FLAGS "-static")
ENDIF (STATIC_BUILD)

ADD_EXECUTABLE(CorrectScatter_multiCore_LMSW CorrectScatter_multiCore_LMSW.cxx)

TARGET_LINK_LIBRARIES(CorrectScatter_multiCore_LMSW ITKCommon
                                    ${ITK_LIBRARIES})

ADD_EXECUTABLE(CorrectAttenuation_SV CorrectAttenuation_SV.cxx)

TARGET_LINK_LIBRARIES(CorrectAttenuation_SV ${VTK_LIBRARIES}
# 				    vtkGraphics vtkRendering vtkHybrid
 				    vtkImaging vtkIO vtkFiltering vtkCommon)

# ItkVtkGlue


ADD_EXECUTABLE(MeshCreator MeshCreator.cxx itkMeshTovtkPolyData.cxx)

TARGET_LINK_LIBRARIES(MeshCreator ITKCommon
                                    vtkGraphics vtkRendering vtkHybrid
                                    vtkImaging vtkIO vtkFiltering vtkCommon

                                    ${ITK_LIBRARIES})

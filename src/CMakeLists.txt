##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of software.
#
# Copyright (c) 2012 University of Pennsylvania. All rights reserved.<br />
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

# ============================================================================
# library target(s)
# ============================================================================

# Add library target for each library using basis_add_library().
#
# This command can not only be used to build libraries from C/C++, but also source
# code written in other programming languages such as Java, Python, Perl,
# MATLAB, and Bash. Note that here we consider modules written in a scripting
# language which are no executables but to be included by other scripts written
# in the particular language as libraries.
#
# Note: Public modules written in a scripting language such as Python, Perl,
#       MATLAB, or Bash which are intended for use by other packages should
#       be placed in the lib/[<lang>/]multiChannelSLIC/ directory,
#       where <lang> is the language name in lowercase and is optional.
#       BASIS will automatically add a library build target for these modules.
basis_add_library (
   supervoxel
   STATIC
   LKM.cpp
   utils.cpp
)

basis_target_link_libraries(supervoxel  ${OpenCV_LIBS}  ${ITK_LIBRARIES})
SET_TARGET_PROPERTIES(supervoxel PROPERTIES COMPILE_FLAGS -fPIC)



# ============================================================================
# executable target(s)
# ============================================================================

# Add executable target for each executable program using basis_add_executable().
#
# This command can not only be used to build executables from C/C++, but also
# source code written in other programming languages such as Java, Python, Perl,
# MATLAB, and Bash.

basis_add_executable (supervoxel_test  supervoxel_test.cpp  LKM.cpp utils.cpp)
basis_target_link_libraries(supervoxel_test ${OpenCV_LIBS}  ${ITK_LIBRARIES})


basis_add_executable (scalarListImgToVecImg   scalarListImgToVecImg.cxx   LKM.cpp utils.cpp)
basis_target_link_libraries(scalarListImgToVecImg   ${ITK_LIBRARIES}  ${OpenCV_LIBS}  )





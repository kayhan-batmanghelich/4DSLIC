##############################################################################
# @file  CMakeLists.txt
# @brief Build configuration of auxiliary data files.
#
# The data within this folder is required to execute the software tools of
# this project. Example data which is only required for an example run of
# the software are, on the other side, available in the example folder.
#
# By default, all content of this directory is simply installed in the
# @c INSTALL_DATA_DIR directory of the installation tree. To group data files,
# place them in proper subdirectories which will be copied to the
# @c DESTINATION directory of the basis_install() function.
#
# Copyright (c) 2014 University of Pennsylvania. All rights reserved.<br />
# See http://www.rad.upenn.edu/sbia/software/license.html or COPYING file.
#
# Contact: SBIA Group <sbia-software at uphs.upenn.edu>
##############################################################################

basis_install_directory ("." "${INSTALL_DATA_DIR}")

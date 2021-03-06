#!/bin/sh
# generate new version header file
exec > new.version.f90
date=`date -R`
host=`hostname -s`
name="$1"
shift
version="$1"
shift
flags="$@"
cat <<EOF
!> Module to provide program version info.
!!
!! This module is autogenerated from a script in order
!! to contain various pieces of information that is only
!! available at compile time. These are the compilation 
!! data and the compiler flags, but also information about
!! the version and the last commit to the source code 
!! management system.
module header

contains

!> Print a program version banner
!! @param channel I/O unit to print banner to.
subroutine version(channel)
  implicit none
  integer, intent(in) :: channel

  write (channel,*) '=================='
  write (channel,*) ' ${name} v${version} '
  write (channel,*) '=================='
  write (channel,*) '-------------------------------------------------------'
  write (channel,*) 'Compile date : ${date} on ${host}'
  write (channel,*) 'Compile flags: ${flags}'
  write (channel,*) '-------------------------------------------------------'
EOF
git log -n 1 --pretty="     write (channel,*) 'Last commit  : %H'"
git log -n 1 --pretty="     write (channel,*) 'Commit date  : %aD'" 
git log -n 1 --pretty="     write (channel,*) 'Commit author: %an <%ae>'" 
cat <<EOF
     write (channel,*) '-------------------------------------------------------'
end subroutine version
end module header
EOF

mv new.version.f90 version.f90

blt_add_executable( NAME serac
                    SOURCES serac.cpp
                    DEPENDS_ON solvers coefficients mfem mpi
                    #DEPENDS_ON solvers coefficients mpi mfem hypre blas lapack ZLIB::ZLIB
                    )

install( TARGETS serac
         RUNTIME DESTINATION bin
         )
